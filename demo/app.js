import init, { scanGrayscale } from "./pkg/zedbar.js";

const dropZone = document.getElementById("drop-zone");
const fileInput = document.getElementById("file-input");
const previewSection = document.getElementById("preview-section");
const previewCanvas = document.getElementById("preview-canvas");
const resultsSection = document.getElementById("results-section");
const resultsContainer = document.getElementById("results");
const resultsSummary = document.getElementById("results-summary");
const noResults = document.getElementById("no-results");
const permalinkSection = document.getElementById("permalink-section");
const permalinkBtn = document.getElementById("permalink-btn");
const errorSection = document.getElementById("error-section");
const errorMessage = document.getElementById("error-message");

let wasmReady = false;
let currentImageUrl = null;

async function initWasm() {
  try {
    await init();
    wasmReady = true;

    // Auto-load image from #url= hash parameter
    const hash = location.hash.slice(1); // remove leading '#'
    const params = new URLSearchParams(hash);
    const url = params.get("url");
    if (url) {
      handleUrl(url);
    }
  } catch (e) {
    showError(`Failed to load WebAssembly module: ${e.message}`);
  }
}

initWasm();

// -- Drop zone events --

dropZone.addEventListener("click", () => fileInput.click());
dropZone.addEventListener("keydown", (e) => {
  if (e.key === "Enter" || e.key === " ") {
    e.preventDefault();
    fileInput.click();
  }
});

fileInput.addEventListener("change", () => {
  if (fileInput.files.length > 0) {
    handleFile(fileInput.files[0]);
  }
});

dropZone.addEventListener("dragover", (e) => {
  e.preventDefault();
  dropZone.classList.add("drag-over");
});

dropZone.addEventListener("dragleave", () => {
  dropZone.classList.remove("drag-over");
});

dropZone.addEventListener("drop", (e) => {
  e.preventDefault();
  dropZone.classList.remove("drag-over");
  const file = e.dataTransfer.files[0];
  if (file) handleFile(file);
});

// Allow paste from clipboard (images or URLs)
document.addEventListener("paste", (e) => {
  const items = e.clipboardData?.items;
  if (!items) return;

  // Prefer a pasted image
  for (const item of items) {
    if (item.type.startsWith("image/")) {
      handleFile(item.getAsFile());
      return;
    }
  }

  // Fall back to pasted text that looks like an image URL
  const text = e.clipboardData.getData("text/plain")?.trim();
  if (text && looksLikeUrl(text)) {
    handleUrl(text);
  }
});

// -- Permalink --

permalinkBtn.addEventListener("click", async () => {
  const permalink = new URL(location.pathname, location.href);
  const hashParams = new URLSearchParams();
  hashParams.set("url", currentImageUrl);
  permalink.hash = hashParams.toString();
  try {
    await navigator.clipboard.writeText(permalink.href);
    permalinkBtn.innerHTML = checkIcon() + " Copied";
    setTimeout(() => {
      permalinkBtn.innerHTML = linkIcon() + " Copy permalink";
    }, 1500);
  } catch {
    // silent fallback — select the URL in the address bar style
  }
});

function fileToDataUrl(file) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(reader.result);
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(file);
  });
}

// -- Main logic --

function looksLikeUrl(text) {
  try {
    const url = new URL(text);
    return url.protocol === "http:" || url.protocol === "https:";
  } catch {
    return false;
  }
}

async function handleUrl(url) {
  if (!wasmReady) {
    showError("WebAssembly module is still loading. Please try again.");
    return;
  }

  hideAll();
  dropZone.classList.add("loading");

  try {
    const resp = await fetch(url);
    if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
    const blob = await resp.blob();
    if (!blob.type.startsWith("image/") && !url.startsWith("data:image/")) {
      throw new Error(`Expected an image, got ${blob.type || "unknown type"}`);
    }
    const bitmap = await createImageBitmap(blob);
    currentImageUrl = url;
    scanBitmap(bitmap);
  } catch (e) {
    showError(`Failed to load image from URL: ${e.message}. Try copying the image itself and pasting it here instead.`);
    dropZone.classList.remove("loading");
  }
}

async function handleFile(file) {
  if (!file.type.startsWith("image/")) {
    showError("Please select an image file.");
    return;
  }

  if (!wasmReady) {
    showError("WebAssembly module is still loading. Please try again.");
    return;
  }

  hideAll();
  dropZone.classList.add("loading");

  try {
    const bitmap = await createImageBitmap(file);
    currentImageUrl = await fileToDataUrl(file);
    scanBitmap(bitmap);
  } catch (e) {
    showError(`Scan failed: ${e}`);
    dropZone.classList.remove("loading");
  }
}

function scanBitmap(bitmap) {
  try {
    const { width, height } = bitmap;

    // Draw to canvas for preview and to extract pixel data
    previewCanvas.width = width;
    previewCanvas.height = height;
    const ctx = previewCanvas.getContext("2d");
    ctx.drawImage(bitmap, 0, 0);
    bitmap.close();

    previewSection.hidden = false;

    // Convert to grayscale
    const imageData = ctx.getImageData(0, 0, width, height);
    const rgba = imageData.data;
    const gray = new Uint8Array(width * height);
    for (let i = 0; i < gray.length; i++) {
      const r = rgba[i * 4];
      const g = rgba[i * 4 + 1];
      const b = rgba[i * 4 + 2];
      // ITU-R BT.601 luma
      gray[i] = (r * 77 + g * 150 + b * 29) >> 8;
    }

    const t0 = performance.now();
    const results = scanGrayscale(gray, width, height);
    const elapsed = performance.now() - t0;
    drawSymbolOverlay(ctx, results);
    displayResults(results, elapsed);
    showPermalink();
  } catch (e) {
    showError(`Scan failed: ${e}`);
  } finally {
    dropZone.classList.remove("loading");
  }
}

function colorForSymbolType(symbolType) {
  switch (symbolType) {
    case "QR-Code":
      return "#22c55e";
    case "SQ-Code":
      return "#3b82f6";
    default:
      return "#f59e0b";
  }
}

// Sort 2D points clockwise around their centroid so a polygon connecting
// them in order forms the convex outline of the symbol. Without this the
// 4 QR corner points (recorded in a non-traversal order) would cross
// themselves when stroked as a polygon.
function pointsAroundCentroid(points) {
  if (points.length <= 2) return points;
  let cx = 0;
  let cy = 0;
  for (const p of points) {
    cx += p.x;
    cy += p.y;
  }
  cx /= points.length;
  cy /= points.length;
  return [...points].sort(
    (a, b) => Math.atan2(a.y - cy, a.x - cx) - Math.atan2(b.y - cy, b.x - cx),
  );
}

// Draws an outline of each detected symbol on top of the source image:
// a polygon for symbols that recorded multiple points (QR corners, scan
// touchpoints) and a bounding rectangle as a fallback. Drawn in image
// coordinates so it scales with the canvas's CSS sizing.
function drawSymbolOverlay(ctx, results) {
  if (!results || results.length === 0) return;
  const stroke = Math.max(2, Math.min(ctx.canvas.width, ctx.canvas.height) / 300);
  ctx.save();
  ctx.lineWidth = stroke;
  ctx.lineJoin = "round";
  ctx.shadowColor = "rgba(0, 0, 0, 0.6)";
  ctx.shadowBlur = stroke * 1.5;
  for (const result of results) {
    ctx.strokeStyle = colorForSymbolType(result.symbolType);
    const points = pointsAroundCentroid(result.points || []);
    if (points.length >= 3) {
      ctx.beginPath();
      ctx.moveTo(points[0].x, points[0].y);
      for (let i = 1; i < points.length; i++) {
        ctx.lineTo(points[i].x, points[i].y);
      }
      ctx.closePath();
      ctx.stroke();
    } else if (result.bounds) {
      const { x, y, width, height } = result.bounds;
      ctx.strokeRect(x, y, width, height);
    }
  }
  ctx.restore();
}

function formatTypeLabel(symbolType, count) {
  // Turn "QR-Code" → "QR code(s)", "CODE-128" → "Code 128(s)", etc.
  const friendly = {
    "QR-Code": "QR code",
    "SQ-Code": "SQ code",
    "CODE-39": "Code 39",
    "CODE-93": "Code 93",
    "CODE-128": "Code 128",
  };
  const label = friendly[symbolType] || symbolType;
  return count === 1 ? label : `${label}s`;
}

function formatElapsed(ms) {
  if (ms < 1) return `${(ms * 1000).toFixed(0)}\u00a0\u00b5s`;
  if (ms < 1000) return `${ms.toFixed(1)}\u00a0ms`;
  return `${(ms / 1000).toFixed(2)}\u00a0s`;
}

function displayResults(results, elapsed) {
  resultsContainer.innerHTML = "";

  if (!results || results.length === 0) {
    noResults.hidden = false;
    resultsSection.hidden = true;
    return;
  }

  resultsSection.hidden = false;

  // Build summary like "2 QR codes, 1 EAN-13"
  const counts = new Map();
  for (const r of results) {
    counts.set(r.symbolType, (counts.get(r.symbolType) || 0) + 1);
  }
  const parts = [];
  for (const [type, count] of counts) {
    const label = formatTypeLabel(type, count);
    parts.push(`${count}\u00a0${label}`);
  }
  const timingStr = `in ${formatElapsed(elapsed)}`;
  const summaryStr = parts.length ? `${parts.join(", ")}` : "";
  resultsSummary.textContent = summaryStr
    ? `\u2014 ${summaryStr} ${timingStr}`
    : `\u2014 ${timingStr}`;

  for (const result of results) {
    const data = result.data;
    const text = result.text;
    const card = document.createElement("div");
    card.className = "result-card";

    const header = document.createElement("div");
    header.className = "result-header";

    const typeBadge = document.createElement("span");
    typeBadge.className = "result-type";
    typeBadge.textContent = result.symbolType;

    const sizeInfo = document.createElement("span");
    sizeInfo.className = "result-size";
    sizeInfo.textContent = `${data.length} byte${data.length !== 1 ? "s" : ""}`;

    header.append(typeBadge, sizeInfo);

    // View tabs
    const tabs = document.createElement("div");
    tabs.className = "view-tabs";

    const hasText = text !== undefined && text !== null;
    const isBinary = !hasText || isBinaryData(data);
    const views = [];

    if (hasText) {
      views.push({ label: "Text", id: "text" });
    }
    views.push({ label: "Hex", id: "hex" });
    views.push({ label: "Base64", id: "base64" });

    for (const view of views) {
      const btn = document.createElement("button");
      btn.className = "view-tab";
      btn.textContent = view.label;
      btn.dataset.view = view.id;
      tabs.appendChild(btn);
    }

    // Data display area
    const dataArea = document.createElement("div");
    dataArea.className = "result-data";

    const content = document.createElement("pre");
    content.className = "data-content";

    const copyBtn = document.createElement("button");
    copyBtn.className = "copy-btn";
    copyBtn.innerHTML = clipboardIcon() + " Copy";

    dataArea.append(content, copyBtn);
    card.append(header, tabs, dataArea);
    resultsContainer.appendChild(card);

    // Tab switching
    function showView(viewId) {
      for (const btn of tabs.querySelectorAll(".view-tab")) {
        btn.classList.toggle("active", btn.dataset.view === viewId);
      }
      content.innerHTML = "";

      if (viewId === "text") {
        content.textContent = text;
      } else if (viewId === "hex") {
        content.appendChild(formatHexDump(data));
      } else if (viewId === "base64") {
        content.textContent = uint8ToBase64(data);
      }

      // Update copy handler
      copyBtn.onclick = () => copyToClipboard(viewId, data, text, copyBtn);
    }

    tabs.addEventListener("click", (e) => {
      const btn = e.target.closest(".view-tab");
      if (btn) showView(btn.dataset.view);
    });

    // Show text tab for text data, hex tab for binary data
    showView(isBinary ? "hex" : views[0].id);
  }
}

// -- Binary detection --

/**
 * Heuristic: data is likely binary if it contains control characters
 * (other than common whitespace: TAB, LF, CR) or null bytes.
 *
 * Intentionally conservative — a QR code containing e.g. a BEL (0x07) in
 * otherwise-text data will be treated as binary, which is acceptable since
 * the text view is still available as a tab.
 */
function isBinaryData(data) {
  for (let i = 0; i < data.length; i++) {
    const b = data[i];
    if (b === 0) return true;
    if (b < 0x20 && b !== 0x09 && b !== 0x0a && b !== 0x0d) return true;
  }
  return false;
}

// -- Hex dump --

function formatHexDump(data) {
  const frag = document.createDocumentFragment();
  const bytesPerLine = 16;

  for (let offset = 0; offset < data.length; offset += bytesPerLine) {
    const line = document.createElement("div");
    line.className = "hex-line";

    const offsetSpan = document.createElement("span");
    offsetSpan.className = "hex-offset";
    offsetSpan.textContent = offset.toString(16).padStart(4, "0") + "  ";

    const bytesSpan = document.createElement("span");
    bytesSpan.className = "hex-bytes";
    const asciiSpan = document.createElement("span");
    asciiSpan.className = "hex-ascii";

    let hexStr = "";
    let asciiStr = "|";

    for (let i = 0; i < bytesPerLine; i++) {
      if (offset + i < data.length) {
        const byte = data[offset + i];
        hexStr += byte.toString(16).padStart(2, "0") + " ";
        asciiStr += byte >= 0x20 && byte <= 0x7e ? String.fromCharCode(byte) : ".";
      } else {
        hexStr += "   ";
        asciiStr += " ";
      }
      if (i === 7) hexStr += " ";
    }

    asciiStr += "|";
    bytesSpan.textContent = hexStr;
    asciiSpan.textContent = asciiStr;

    line.append(offsetSpan, bytesSpan, asciiSpan);
    frag.appendChild(line);
  }

  return frag;
}

// -- Copy --

async function copyToClipboard(viewId, data, text, btn) {
  let copyText;

  if (viewId === "text") {
    copyText = text;
  } else if (viewId === "hex") {
    copyText = formatHexString(data);
  } else if (viewId === "base64") {
    copyText = uint8ToBase64(data);
  }

  try {
    await navigator.clipboard.writeText(copyText);
    btn.innerHTML = checkIcon() + " Copied";
    btn.classList.add("copied");
    setTimeout(() => {
      btn.innerHTML = clipboardIcon() + " Copy";
      btn.classList.remove("copied");
    }, 1500);
  } catch {
    // Fallback
    const ta = document.createElement("textarea");
    ta.value = copyText;
    ta.style.position = "fixed";
    ta.style.opacity = "0";
    document.body.appendChild(ta);
    ta.select();
    document.execCommand("copy");
    document.body.removeChild(ta);
    btn.innerHTML = checkIcon() + " Copied";
    btn.classList.add("copied");
    setTimeout(() => {
      btn.innerHTML = clipboardIcon() + " Copy";
      btn.classList.remove("copied");
    }, 1500);
  }
}

function formatHexString(data) {
  const bytesPerLine = 16;
  const lines = [];
  for (let offset = 0; offset < data.length; offset += bytesPerLine) {
    let hexStr = offset.toString(16).padStart(4, "0") + "  ";
    let asciiStr = "|";
    for (let i = 0; i < bytesPerLine; i++) {
      if (offset + i < data.length) {
        const byte = data[offset + i];
        hexStr += byte.toString(16).padStart(2, "0") + " ";
        asciiStr += byte >= 0x20 && byte <= 0x7e ? String.fromCharCode(byte) : ".";
      } else {
        hexStr += "   ";
        asciiStr += " ";
      }
      if (i === 7) hexStr += " ";
    }
    asciiStr += "|";
    lines.push(hexStr + asciiStr);
  }
  return lines.join("\n");
}

function uint8ToBase64(data) {
  let binary = "";
  for (const byte of data) {
    binary += String.fromCharCode(byte);
  }
  return btoa(binary);
}

// -- Utility --

function showPermalink() {
  if (currentImageUrl) {
    permalinkBtn.innerHTML = linkIcon() + " Copy permalink";
    permalinkSection.hidden = false;
  }
}

function hideAll() {
  resultsSection.hidden = true;
  noResults.hidden = true;
  errorSection.hidden = true;
  permalinkSection.hidden = true;
}

function showError(msg) {
  hideAll();
  errorMessage.textContent = msg;
  errorSection.hidden = false;
}

function clipboardIcon() {
  return `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>`;
}

function linkIcon() {
  return `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M10 13a5 5 0 0 0 7.54.54l3-3a5 5 0 0 0-7.07-7.07l-1.72 1.71"/><path d="M14 11a5 5 0 0 0-7.54-.54l-3 3a5 5 0 0 0 7.07 7.07l1.71-1.71"/></svg>`;
}

function checkIcon() {
  return `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="20 6 9 17 4 12"/></svg>`;
}

// -- Install tabs --

const installTabs = document.querySelector(".install-tabs");
if (installTabs) {
  installTabs.addEventListener("click", (e) => {
    const tab = e.target.closest(".install-tab");
    if (!tab) return;

    const lang = tab.dataset.lang;

    // Update active tab
    for (const t of installTabs.querySelectorAll(".install-tab")) {
      t.classList.toggle("active", t.dataset.lang === lang);
    }

    // Show/hide content
    for (const content of document.querySelectorAll(".install-content")) {
      content.hidden = content.dataset.lang !== lang;
    }
  });
}
