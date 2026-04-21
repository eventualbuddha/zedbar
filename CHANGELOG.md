# Changelog

## [0.2.5](https://github.com/eventualbuddha/zedbar/compare/v0.2.4...v0.2.5) (2026-04-21)


### Bug Fixes

* detect small QR codes in heavily noisy images via multi-region retry ([#45](https://github.com/eventualbuddha/zedbar/issues/45)) ([898b671](https://github.com/eventualbuddha/zedbar/commit/898b67124fc9be2871a54ad2c6d18c315664d881))
* enable `retry_undecoded_regions` in zedbarimg for small QR codes ([#43](https://github.com/eventualbuddha/zedbar/issues/43)) ([0319202](https://github.com/eventualbuddha/zedbar/commit/03192020f146150f508ea6d87f4f2bb0d6c447f0))

## [0.2.4](https://github.com/eventualbuddha/zedbar/compare/v0.2.3...v0.2.4) (2026-04-19)


### Bug Fixes

* correct EAN-8 digit reversal and QR finder buffer corruption ([#40](https://github.com/eventualbuddha/zedbar/issues/40)) ([e34ef7f](https://github.com/eventualbuddha/zedbar/commit/e34ef7fb3868684eafa06e8d6cf1e856314196d6))

## [0.2.3](https://github.com/eventualbuddha/zedbar/compare/v0.2.2...v0.2.3) (2026-04-18)


### Bug Fixes

* allow reading binary QR/SQ code data ([#37](https://github.com/eventualbuddha/zedbar/issues/37)) ([ca5ca98](https://github.com/eventualbuddha/zedbar/commit/ca5ca9840c4a6f722eb7bf39dc6524e6afd2a6bb))

## [0.2.2](https://github.com/eventualbuddha/zedbar/compare/v0.2.1...v0.2.2) (2026-04-18)


### Bug Fixes

* detect small QR codes in large images via crop+upscale ([#30](https://github.com/eventualbuddha/zedbar/issues/30)) ([c545cb2](https://github.com/eventualbuddha/zedbar/commit/c545cb264dd19f9521e3fa49cdb1844a0e0ba12c))
* improve QR detection for images with many QR codes ([#33](https://github.com/eventualbuddha/zedbar/issues/33)) ([b503f35](https://github.com/eventualbuddha/zedbar/commit/b503f3528f42716926b01c0261b4e1983cfadb79))

## [0.2.1](https://github.com/eventualbuddha/zedbar/compare/v0.2.0...v0.2.1) (2026-02-11)


### Bug Fixes

* remove deprecated command parameter from release-please v4 ([#16](https://github.com/eventualbuddha/zedbar/issues/16)) ([6532d1e](https://github.com/eventualbuddha/zedbar/commit/6532d1e7467746cc43dc494243fead25d373783f))
* simplify release-please to use release-type: rust ([#17](https://github.com/eventualbuddha/zedbar/issues/17)) ([cb11a9d](https://github.com/eventualbuddha/zedbar/commit/cb11a9d3102225efe8ca9ad9e4a59545e96c9cac))

## [0.2.0](https://github.com/eventualbuddha/zedbar/compare/v0.1.1...v0.2.0) (2026-02-07)


### Features

* add zedbarimg CLI for npm with WebP support ([#13](https://github.com/eventualbuddha/zedbar/issues/13)) ([5c034dc](https://github.com/eventualbuddha/zedbar/commit/5c034dcb313ae0412399270ffcd9342941f696a3))


### Bug Fixes

* remove deprecated command parameter from release-please v4 ([#16](https://github.com/eventualbuddha/zedbar/issues/16)) ([6532d1e](https://github.com/eventualbuddha/zedbar/commit/6532d1e7467746cc43dc494243fead25d373783f))
* simplify release-please to use release-type: rust ([#17](https://github.com/eventualbuddha/zedbar/issues/17)) ([cb11a9d](https://github.com/eventualbuddha/zedbar/commit/cb11a9d3102225efe8ca9ad9e4a59545e96c9cac))
