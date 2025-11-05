# Implementation Guide: Removing the `userdata` Pointer

This document provides a step-by-step guide for implementing the recommended solution (return values instead of callbacks) to eliminate the unsafe `userdata` pointer.

## Quick Reference

**Current State**: Decoder calls back to ImageScanner via raw pointer
**Target State**: Decoder returns symbol data, ImageScanner handles it
**Files to modify**: 3 main files
**Estimated complexity**: Medium (requires careful refactoring but straightforward)

## Step-by-Step Implementation

### Step 1: Define Symbol Data Structure

Create a structure to hold decoded symbol information that will be returned.

**File**: `src/decoder.rs`
**Add after line 632** (before `zbar_decoder_t` definition):

```rust
/// Data for a decoded symbol ready to be added to results
#[derive(Debug)]
pub(crate) struct DecodedSymbolData {
    pub(crate) symbol_type: SymbolType,
    pub(crate) data: Vec<u8>,
    pub(crate) modifiers: c_uint,
    pub(crate) direction: c_int,
}
```

### Step 2: Add Pending Symbol to Decoder

Add a field to store the decoded symbol data before returning it.

**File**: `src/decoder.rs`
**Modify** `zbar_decoder_t` struct (around line 635):

```rust
pub(crate) struct zbar_decoder_t {
    // ... existing fields ...

    buffer: Vec<u8>,
-   pub(crate) userdata: *mut c_void,  // ← REMOVE THIS LINE

+   // Pending symbol data to return
+   pending_symbol: Option<DecodedSymbolData>,

    // Configuration (new type-safe system)
    pub(crate) config: DecoderState,
    // ... rest of fields ...
}
```

**Update Default implementation** (around line 688):

```rust
fn with_config(config: DecoderState) -> Self {
    let mut decoder = Self {
        // ... existing initialization ...
        buffer: Vec::with_capacity(BUFFER_MIN as usize),
-       userdata: null_mut(),  // ← REMOVE THIS LINE
+       pending_symbol: None,  // ← ADD THIS LINE
        config: config.clone(),
        // ... rest of initialization ...
    };
    // ...
}
```

### Step 3: Change decode_width to Return Symbol Data

Modify `decode_width` to return the pending symbol instead of calling the handler.

**File**: `src/decoder.rs`
**Modify** the `decode_width` method (line 926):

```rust
-pub(crate) unsafe fn decode_width(&mut self, w: c_uint) -> SymbolType {
+pub(crate) fn decode_width(&mut self, w: c_uint) -> Option<DecodedSymbolData> {
+    // Clear any previous pending symbol
+    self.pending_symbol = None;

     let mut sym = SymbolType::None;

     // Store width in circular buffer
     self.w[(self.idx & (DECODE_WINDOW - 1) as u8) as usize] = w;

     // ... existing decoding logic ...

     self.idx = self.idx.wrapping_add(1);
     self.type_ = sym;

     if sym != SymbolType::None {
         if self.lock != SymbolType::None
             && sym > SymbolType::Partial
             && sym != SymbolType::QrCode
         {
             self.release_lock(sym);
         }

-        symbol_handler(self);  // ← REMOVE THIS LINE
+        // Instead of calling handler, prepare symbol data to return
+        if sym > SymbolType::Partial && sym != SymbolType::QrCode {
+            self.pending_symbol = Some(DecodedSymbolData {
+                symbol_type: sym,
+                data: self.buffer_slice().to_vec(),
+                modifiers: self.modifiers,
+                direction: self.direction,
+            });
+        }
     }

-    sym
+    self.pending_symbol.take()
 }
```

### Step 4: Remove Userdata Methods

Remove the getter and setter for userdata.

**File**: `src/decoder.rs`
**Delete** lines 1014-1022:

```rust
-    /// Get user data pointer
-    pub(crate) fn get_userdata(&self) -> *mut c_void {
-        self.userdata
-    }
-
-    /// Set user data pointer
-    pub(crate) fn set_userdata(&mut self, userdata: *mut c_void) {
-        self.userdata = userdata;
-    }
```

### Step 5: Update LineScanner to Propagate Symbols

Modify the line scanner to pass through symbol data.

**File**: `src/line_scanner.rs`
**Modify** `process_edge` (line 253):

```rust
-fn process_edge(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
+fn process_edge(&mut self, decoder: &mut zbar_decoder_t) -> Option<DecodedSymbolData> {
     // ... existing edge calculation ...

     self.width = self.cur_edge - self.last_edge;
     self.last_edge = self.cur_edge;

     let width = self.width;
-    unsafe { decoder.decode_width(width) }
+    decoder.decode_width(width)
 }
```

**Modify** `scan_y` (line 144):

```rust
-pub(crate) fn scan_y(&mut self, y: c_int, decoder: &mut zbar_decoder_t) -> SymbolType {
+pub(crate) fn scan_y(&mut self, y: c_int, decoder: &mut zbar_decoder_t)
+    -> Option<DecodedSymbolData>
+{
     // ... existing pixel processing ...

-    let mut edge = SymbolType::None;
+    let mut edge = None;

     // 2nd zero-crossing is 1st local min/max - could be edge
     if (y2_1 == 0 || ((y2_1 > 0) == (y2_2 < 0))) && (self.calc_thresh() <= y1_1.unsigned_abs())
     {
         // ... existing logic ...

         if y1_rev {
             // intensity change reversal - finalize previous edge
             edge = self.process_edge(decoder);
         }

         // ... rest of logic ...
     }

     self.x = x + 1;
     edge
 }
```

**Modify** `scanner_flush` (line 99):

```rust
-pub(crate) fn scanner_flush(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
+pub(crate) fn scanner_flush(&mut self, decoder: &mut zbar_decoder_t)
+    -> Option<DecodedSymbolData>
+{
     if self.y1_sign == 0 {
-        return SymbolType::None;
+        return None;
     }

     // ... existing logic ...

     self.y1_sign = 0;
     self.width = 0;
-    unsafe { decoder.decode_width(0) }
+    decoder.decode_width(0)
 }
```

**Modify** `new_scan` (line 120):

```rust
-pub(crate) fn new_scan(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
-    let mut edge = SymbolType::None;
+pub(crate) fn new_scan(&mut self, decoder: &mut zbar_decoder_t) -> Vec<DecodedSymbolData> {
+    let mut symbols = Vec::new();

     while self.y1_sign != 0 {
-        let tmp = self.scanner_flush(decoder);
-        if tmp > edge {
-            edge = tmp;
+        if let Some(symbol) = self.scanner_flush(decoder) {
+            symbols.push(symbol);
         }
     }

     // reset scanner and associated decoder
     // ... existing reset logic ...

     decoder.new_scan();
-    edge
+    symbols
 }
```

**Modify** `quiet_border` (line 224):

```rust
-pub(crate) fn quiet_border(&mut self, decoder: &mut zbar_decoder_t) {
+pub(crate) fn quiet_border(&mut self, decoder: &mut zbar_decoder_t)
+    -> Vec<DecodedSymbolData>
+{
     // Flush scanner pipeline twice
     self.scanner_flush(decoder);
     self.scanner_flush(decoder);

     // Start new scan
-    self.new_scan(decoder);
+    self.new_scan(decoder)
 }
```

### Step 6: Update ImageScanner Scan Loop

Update the image scanner to handle returned symbols.

**File**: `src/img_scanner.rs`
**Remove** the userdata setup (delete lines 232-234):

```rust
-    // Set up decoder's back-pointer to this scanner for symbol callbacks
-    let scanner_ptr = self as *mut _ as *mut c_void;
-    self.dcode.set_userdata(scanner_ptr);
```

**Update** horizontal scan loop (around line 274-280):

```rust
 while (x as c_uint) < w {
     let d = data[p as usize];
     x += 1;
     p += 1;
-    self.scn.scan_y(d as c_int, &mut self.dcode);
+    if let Some(symbol_data) = self.scn.scan_y(d as c_int, &mut self.dcode) {
+        self.handle_decoded_symbol(symbol_data);
+    }
 }
-self.scn.quiet_border(&mut self.dcode);
+for symbol_data in self.scn.quiet_border(&mut self.dcode) {
+    self.handle_decoded_symbol(symbol_data);
+}
```

**Update** the second horizontal pass (around line 294-300):

```rust
 while x >= 0 {
     let d = data[p as usize];
     x -= 1;
     p -= 1;
-    self.scn.scan_y(d as c_int, &mut self.dcode);
+    if let Some(symbol_data) = self.scn.scan_y(d as c_int, &mut self.dcode) {
+        self.handle_decoded_symbol(symbol_data);
+    }
 }
-self.scn.quiet_border(&mut self.dcode);
+for symbol_data in self.scn.quiet_border(&mut self.dcode) {
+    self.handle_decoded_symbol(symbol_data);
+}
```

**Update** vertical scan loops similarly (around lines 332-365).

### Step 7: Add Symbol Handler as Method

Convert the standalone `symbol_handler` function into a method.

**File**: `src/img_scanner.rs`
**Replace** `symbol_handler` function (delete lines 472-538) with a method in the `impl` block:

```rust
impl zbar_image_scanner_t {
    // ... existing methods ...

+    /// Handle a decoded symbol by adding it to results or updating duplicates
+    fn handle_decoded_symbol(&mut self, symbol_data: DecodedSymbolData) {
+        let symbol_type = symbol_data.symbol_type;
+
+        // Calculate position if position tracking is enabled
+        let (x, y) = if self.config.position_tracking {
+            let scn = &self.scn;
+            let w = scn.width();
+            let u = self.umin + self.du * scn.get_edge(w, 0) as c_int;
+            if self.dx != 0 {
+                (u, self.v)
+            } else {
+                (self.v, u)
+            }
+        } else {
+            (0, 0)
+        };
+
+        // Check for duplicate
+        if let Some(sym) = self.find_duplicate_symbol(symbol_type, &symbol_data.data) {
+            sym.quality += 1;
+            if self.config.position_tracking {
+                sym.add_point(x, y);
+            }
+            return;
+        }
+
+        // Create new symbol
+        let mut sym = zbar_symbol_t::new(symbol_type);
+        sym.modifiers = symbol_data.modifiers;
+        sym.data.extend_from_slice(&symbol_data.data);
+
+        if self.config.position_tracking {
+            sym.add_point(x, y);
+        }
+
+        // Set orientation
+        let dir = symbol_data.direction;
+        if dir != 0 {
+            sym.orient = (if self.dy != 0 { 1 } else { 0 }) + ((self.du ^ dir) & 2);
+        }
+
+        self.add_symbol(sym);
+    }
}
```

**Delete** the old `symbol_handler` function entirely.

### Step 8: Update Imports

**File**: `src/decoder.rs` (top of file):

Remove the import of `symbol_handler`:
```rust
-use crate::img_scanner::symbol_handler;
```

Export the new type:
```rust
+pub(crate) use decoder::DecodedSymbolData;
```

**File**: `src/img_scanner.rs`:

Remove now-unused imports if any.

### Step 9: Testing

After making all changes:

1. **Compile**:
   ```bash
   cargo build
   ```

2. **Run tests**:
   ```bash
   cargo test
   ```

3. **Check for unsafe**:
   ```bash
   rg "unsafe" src/
   ```

   Should only show unsafe in:
   - FFI boundary functions
   - Specific performance-critical sections with clear safety comments
   - NOT in symbol_handler (which should be deleted)

4. **Run clippy**:
   ```bash
   cargo clippy
   ```

## Checklist

- [ ] Created `DecodedSymbolData` struct
- [ ] Removed `userdata` field from `zbar_decoder_t`
- [ ] Added `pending_symbol` field to `zbar_decoder_t`
- [ ] Changed `decode_width` to return `Option<DecodedSymbolData>`
- [ ] Removed `symbol_handler` call from `decode_width`
- [ ] Removed `get_userdata()` and `set_userdata()` methods
- [ ] Updated `process_edge` in `line_scanner.rs`
- [ ] Updated `scan_y` in `line_scanner.rs`
- [ ] Updated `scanner_flush` in `line_scanner.rs`
- [ ] Updated `new_scan` in `line_scanner.rs`
- [ ] Updated `quiet_border` in `line_scanner.rs`
- [ ] Removed `set_userdata()` call in `scan_image_internal`
- [ ] Updated all `scan_y` calls in horizontal scan loops
- [ ] Updated all `scan_y` calls in vertical scan loops
- [ ] Updated all `quiet_border` calls
- [ ] Converted `symbol_handler` to `handle_decoded_symbol` method
- [ ] Updated imports
- [ ] Code compiles without errors
- [ ] All tests pass
- [ ] No unsafe code in symbol handling
- [ ] Clippy warnings addressed

## Expected Benefits

After implementation:

1. ✅ **No unsafe code** for symbol handling
2. ✅ **No raw pointers** creating circular references
3. ✅ **Clear data flow** - symbols flow up the call stack
4. ✅ **Better testability** - can test decoder in isolation
5. ✅ **More maintainable** - explicit rather than callback-based
6. ✅ **Rust idiomatic** - follows ownership principles

## Troubleshooting

### "Cannot return value referencing local variable"

Make sure `DecodedSymbolData` contains **owned** data (`Vec<u8>`), not references.

### "Borrow checker errors in scan loops"

Ensure you're not holding references across the `handle_decoded_symbol` call. The symbol data should be moved/consumed.

### "Tests failing"

Check that QR code handling still works - QR codes are handled separately and shouldn't go through the new symbol data path (they return early in the decoder).

## Alternative: Smaller First Step

If the full refactor is too large, you can start with just the 1D decoders:

1. Keep QR handling as-is (returns `SymbolType::QrCode`)
2. Only return `DecodedSymbolData` for 1D symbologies
3. Make `decode_width` return `Result<Option<DecodedSymbolData>, SymbolType>`
4. Handle `Ok(Some(data))` as new path, `Err(QrCode)` as old path

This allows incremental migration and testing.
