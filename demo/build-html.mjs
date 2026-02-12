#!/usr/bin/env node
/**
 * Build script that processes index.src.html and generates index.html
 * with syntax-highlighted code blocks using shiki.
 */

import { readFile, writeFile } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';
import { createHighlighter } from 'shiki';

const __dirname = dirname(fileURLToPath(import.meta.url));
const projectDir = join(__dirname, '..');

async function main() {
  // Read version from Cargo.toml
  const cargoToml = await readFile(join(projectDir, 'Cargo.toml'), 'utf-8');
  const versionMatch = cargoToml.match(/^version\s*=\s*"([^"]+)"/m);
  if (!versionMatch) {
    throw new Error('Could not find version in Cargo.toml');
  }
  const fullVersion = versionMatch[1];
  // Use major.minor for the dependency specifier (e.g. "0.2.1" -> "0.2")
  const depVersion = fullVersion.split('.').slice(0, 2).join('.');

  const highlighter = await createHighlighter({
    themes: ['github-light', 'one-dark-pro'],
    langs: ['rust', 'javascript', 'bash', 'toml'],
  });

  let html = await readFile(join(__dirname, 'index.src.html'), 'utf-8');

  // Replace version placeholder
  html = html.replace(/\{\{ZEDBAR_VERSION\}\}/g, depVersion);

  // Find all <code data-lang="...">...</code> blocks and highlight them
  const codeBlockRegex = /<code data-lang="(\w+)">([\s\S]*?)<\/code>/g;

  html = html.replace(codeBlockRegex, (match, lang, code) => {
    // Decode HTML entities back to actual characters for highlighting
    const decoded = code
      .replace(/&amp;/g, '&')
      .replace(/&lt;/g, '<')
      .replace(/&gt;/g, '>')
      .replace(/&quot;/g, '"');

    const highlighted = highlighter.codeToHtml(decoded, {
      lang,
      themes: {
        light: 'github-light',
        dark: 'one-dark-pro',
      },
    });

    // Extract just the inner content from shiki's output
    // shiki wraps in <pre class="shiki"><code>...</code></pre>
    // We want to keep our own <pre> and just replace the <code> content
    const innerMatch = highlighted.match(/<code>([\s\S]*)<\/code>/);
    if (innerMatch) {
      return `<code>${innerMatch[1]}</code>`;
    }
    return match;
  });

  // CSS for dark mode - light mode uses inline color: styles directly
  const css = `
/* Shiki syntax highlighting - dark mode override */
@media (prefers-color-scheme: dark) {
  pre code span[style] { color: var(--shiki-dark) !important; }
}`;

  // Inject CSS before </head>
  html = html.replace('</head>', `  <style>${css}\n  </style>\n</head>`);

  await writeFile(join(__dirname, 'index.html'), html);
  console.log('Generated index.html with syntax highlighting');
}

await main();
