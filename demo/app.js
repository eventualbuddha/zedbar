import init, { scan_grayscale } from "./pkg/zedbar.js";

const dropZone = document.getElementById("drop-zone");
const fileInput = document.getElementById("file-input");
const previewSection = document.getElementById("preview-section");
const previewCanvas = document.getElementById("preview-canvas");
const resultsSection = document.getElementById("results-section");
const resultsContainer = document.getElementById("results");
const noResults = document.getElementById("no-results");
const errorSection = document.getElementById("error-section");
const errorMessage = document.getElementById("error-message");

let wasmReady = false;

async function initWasm() {
  try {
    await init();
    wasmReady = true;
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

// Allow paste from clipboard
document.addEventListener("paste", (e) => {
  const items = e.clipboardData?.items;
  if (!items) return;
  for (const item of items) {
    if (item.type.startsWith("image/")) {
      handleFile(item.getAsFile());
      return;
    }
  }
});

// -- Main logic --

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

    const results = scan_grayscale(gray, width, height);
    displayResults(results);
  } catch (e) {
    showError(`Scan failed: ${e}`);
  } finally {
    dropZone.classList.remove("loading");
  }
}

function displayResults(results) {
  resultsContainer.innerHTML = "";

  if (!results || results.length === 0) {
    noResults.hidden = false;
    resultsSection.hidden = true;
    return;
  }

  resultsSection.hidden = false;

  for (const result of results) {
    const data = result.data;
    const text = result.text;
    const card = document.createElement("div");
    card.className = "result-card";

    const header = document.createElement("div");
    header.className = "result-header";

    const typeBadge = document.createElement("span");
    typeBadge.className = "result-type";
    typeBadge.textContent = result.symbol_type;

    const sizeInfo = document.createElement("span");
    sizeInfo.className = "result-size";
    sizeInfo.textContent = `${data.length} byte${data.length !== 1 ? "s" : ""}`;

    header.append(typeBadge, sizeInfo);

    // View tabs
    const tabs = document.createElement("div");
    tabs.className = "view-tabs";

    const views = [];

    if (text !== undefined && text !== null) {
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

    // Show first tab
    showView(views[0].id);
  }
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

function hideAll() {
  resultsSection.hidden = true;
  noResults.hidden = true;
  errorSection.hidden = true;
}

function showError(msg) {
  hideAll();
  errorMessage.textContent = msg;
  errorSection.hidden = false;
}

function clipboardIcon() {
  return `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/></svg>`;
}

function checkIcon() {
  return `<svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polyline points="20 6 9 17 4 12"/></svg>`;
}
