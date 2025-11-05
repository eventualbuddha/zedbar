# Symbol Decoding Data Flow

This document provides detailed diagrams showing how symbols are decoded and the specific path of the unsafe `userdata` pointer.

## Object Layout and Ownership

```mermaid
classDiagram
    class ImageScanner {
        -LineScanner scn
        -Decoder dcode
        -QrReader qr
        -SymbolSet syms
        -ImageScannerConfig config
        +scan_image(img) Result~i32~
        +add_symbol(sym)
        +find_duplicate_symbol() Option
        -scan_image_internal()
    }

    class LineScanner {
        -uint x
        -uint width
        -int[4] y0
        +scan_y(pixel, decoder) SymbolType
        +process_edge(decoder) SymbolType
        +new_scan(decoder)
    }

    class Decoder {
        -uint[16] w
        -SymbolType type_
        -Vec~u8~ buffer
        -void* userdata ⚠️
        -EanDecoder ean
        -I25Decoder i25
        -Code39Decoder code39
        +decode_width(width) SymbolType
        +get_userdata() void*
        +set_userdata(ptr)
    }

    class SymbolSet {
        +Vec~Symbol~ symbols
    }

    ImageScanner "1" --> "1" LineScanner: owns
    ImageScanner "1" --> "1" Decoder: owns
    ImageScanner "1" --> "1" SymbolSet: owns
    Decoder "1" ..> "1" ImageScanner: userdata (UNSAFE!)

    note for Decoder "The userdata pointer creates a circular reference that violates Rust's borrow rules"
```

## Complete Call Chain with Borrow States

```mermaid
sequenceDiagram
    autonumber
    participant App
    participant IS as ImageScanner
    participant LS as LineScanner
    participant Dec as Decoder
    participant SH as symbol_handler

    Note over IS: State: Not borrowed
    App->>+IS: scan_image(&mut self, img)
    Note over IS: State: &mut borrow (ACTIVE)

    IS->>IS: scan_image_internal()
    Note over IS,Dec: Setup phase
    IS->>Dec: set_userdata(self as *mut void)
    Note over Dec: userdata now points to IS

    IS->>LS: new_scan(&mut self, &mut decoder)

    loop Each horizontal scan line
        loop Each pixel in line
            Note over IS: Still &mut borrowed!
            IS->>+LS: scan_y(pixel, &mut self.dcode)
            Note over LS: &mut borrow of decoder

            LS->>LS: Detect edge?
            alt Edge detected
                LS->>+Dec: decode_width(&mut self, width)
                Note over Dec: &mut borrow of decoder

                Dec->>Dec: Try decode EAN
                Dec->>Dec: Try decode Code39
                Dec->>Dec: Try decode Code128
                Note over Dec: etc...

                alt Symbol decoded successfully
                    Note over Dec,SH: THE PROBLEM STARTS HERE
                    Dec->>+SH: symbol_handler(&mut decoder)

                    SH->>Dec: get_userdata() -> *mut void
                    Note over SH: Cast to *mut ImageScanner

                    SH->>IS: (*iscn).find_duplicate_symbol()
                    Note over IS: ⚠️ VIOLATION!<br/>Already &mut borrowed<br/>Creating 2nd mutable reference

                    alt Duplicate found
                        SH->>IS: (*iscn).syms[i].quality += 1
                        SH->>IS: (*iscn).syms[i].add_point()
                    else New symbol
                        SH->>IS: (*iscn).add_symbol(new_sym)
                        Note over IS: Pushes to syms.symbols Vec
                    end

                    SH-->>-Dec: return
                end

                Dec-->>-LS: return SymbolType
            end
            LS-->>-IS: return
        end
        IS->>LS: quiet_border(&mut decoder)
    end

    IS-->>-App: Ok(num_symbols)
    Note over IS: State: Not borrowed
```

## Memory Layout Showing the Circular Reference

```mermaid
graph LR
    subgraph "Stack Frame: scan_image"
        SF[" &mut self "]
    end

    subgraph "Heap: ImageScanner @ 0x1000"
        IS_Header["ImageScanner"]
        IS_scn["scn: LineScanner"]
        IS_dcode["dcode: Decoder"]
        IS_syms["syms: SymbolSet"]
    end

    subgraph "Inside Decoder @ 0x1000 + offset"
        Dec_w["w: [u32; 16]"]
        Dec_buffer["buffer: Vec&lt;u8&gt;"]
        Dec_userdata["userdata: *mut void<br/>= 0x1000"]
        Dec_ean["ean: EanDecoder"]
    end

    subgraph "Inside SymbolSet"
        Syms_vec["symbols: Vec&lt;Symbol&gt;"]
    end

    SF -.->|"Mutable borrow"| IS_Header
    IS_Header --> IS_scn
    IS_Header --> IS_dcode
    IS_Header --> IS_syms

    IS_dcode --> Dec_w
    IS_dcode --> Dec_buffer
    IS_dcode --> Dec_userdata
    IS_dcode --> Dec_ean

    IS_syms --> Syms_vec

    Dec_userdata -.->|"Raw pointer<br/>back to 0x1000"| IS_Header

    style SF fill:#e1f5ff
    style Dec_userdata fill:#ffcccc
    style IS_Header fill:#ffcccc
    style Syms_vec fill:#ffcccc

    Note1["⚠️ PROBLEM:<br/>scan_image has &mut ImageScanner<br/>symbol_handler creates 2nd &mut via raw pointer<br/>Both can access syms simultaneously"]
    style Note1 fill:#ffe6e6
```

## Specific Code Locations

### 1. Setup (img_scanner.rs:232-234)

```rust
fn scan_image_internal(&mut self, img: &mut zbar_image_t) -> Option<zbar_symbol_set_t> {
    // Set up decoder's back-pointer to this scanner for symbol callbacks
    let scanner_ptr = self as *mut _ as *mut c_void;  // ← Creates raw pointer
    self.dcode.set_userdata(scanner_ptr);              // ← Stores in decoder
    // ...
}
```

At this point:
- `self` is borrowed mutably for the entire function
- `scanner_ptr` is a raw pointer to the same memory
- Both exist simultaneously!

### 2. The Scan Loop (img_scanner.rs:274-279)

```rust
while (x as c_uint) < w {
    let d = data[p as usize];
    x += 1;
    p += 1;
    self.scn.scan_y(d as c_int, &mut self.dcode);  // ← Passes &mut decoder
}
```

Still inside the `&mut self` borrow from `scan_image_internal`.

### 3. Edge Detection (line_scanner.rs:253-266)

```rust
fn process_edge(&mut self, decoder: &mut zbar_decoder_t) -> SymbolType {
    // ... edge calculation ...
    self.width = self.cur_edge - self.last_edge;
    self.last_edge = self.cur_edge;

    let width = self.width;
    unsafe { decoder.decode_width(width) }  // ← Calls into decoder
}
```

### 4. Decoding (decoder.rs:926-1008)

```rust
pub(crate) unsafe fn decode_width(&mut self, w: c_uint) -> SymbolType {
    // ... store width ...

    // Try each decoder
    if self.ean.enable {
        let tmp = zbar_decode_ean(&mut *self);
        if tmp != SymbolType::None {
            sym = tmp;
        }
    }
    // ... more decoders ...

    if sym != SymbolType::None {
        // ...
        symbol_handler(self);  // ← THE CALLBACK
    }

    sym
}
```

### 5. The Unsafe Callback (img_scanner.rs:479-538)

```rust
pub(crate) unsafe fn symbol_handler(dcode: &mut zbar_decoder_t) {
    let iscn = dcode.get_userdata() as *mut zbar_image_scanner_t;  // ← Get raw pointer
    // ...

    let data = dcode.buffer_slice();
    if let Some(sym) = (*iscn).find_duplicate_symbol(symbol_type, data) {  // ← Deref!
        sym.quality += 1;  // ← Mutate through raw pointer
        if (*iscn).config.position_tracking {
            sym.add_point(x, y);
        }
        return;
    }

    // ... or create new symbol ...
    (*iscn).add_symbol(sym);  // ← Mutate through raw pointer
}
```

Here's the problem:
- `(*iscn)` creates a mutable reference from the raw pointer
- But `iscn` points to the same ImageScanner that's already borrowed by `scan_image_internal`
- We now have two mutable references to the same memory!

## Why This is Undefined Behavior

According to Rust's rules:

1. **Aliasing XOR Mutability**: You can have either:
   - Many immutable references, OR
   - Exactly ONE mutable reference

2. **Our violation**:
   - `scan_image_internal(&mut self)` creates mutable reference #1
   - `(*iscn)` in `symbol_handler` creates mutable reference #2
   - Both are active simultaneously
   - This is undefined behavior, even in unsafe code

3. **Why it "works"**:
   - The compiler can't prove the aliasing due to the raw pointer
   - No actual memory corruption happens in practice (yet)
   - But the optimizer could theoretically break this code

## Visual: The Borrow Stack Violation

```
┌─────────────────────────────────────────┐
│ Active Borrows (what Rust sees)         │
├─────────────────────────────────────────┤
│ &mut ImageScanner (scan_image)          │  ← Borrow 1
│   └─ &mut Decoder (scan_y)              │
│       └─ decode_width running...        │
├─────────────────────────────────────────┤
│ Raw pointer manipulation (hidden)       │
├─────────────────────────────────────────┤
│ &mut ImageScanner (from userdata)       │  ← Borrow 2 (ILLEGAL!)
│   └─ Calling add_symbol()               │
│   └─ Mutating syms.symbols              │
└─────────────────────────────────────────┘
      ↑                    ↑
      └────── SAME MEMORY ─┘

This violates Rust's aliasing rules!
```

## Summary

The `userdata` pointer creates a hidden circular reference that bypasses Rust's borrow checker:

1. **Outer context**: `ImageScanner` has a mutable borrow active
2. **Stored pointer**: Decoder stores a raw pointer to `ImageScanner`
3. **Callback**: Decoder calls back through the pointer
4. **Violation**: Creates a second mutable reference while the first is still active

The only way to fix this properly is to eliminate the callback pattern and use one of the solutions outlined in `UNSAFE_USERDATA_PROBLEM.md`.
