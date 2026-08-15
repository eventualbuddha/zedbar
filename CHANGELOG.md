# Changelog

## [0.5.1](https://github.com/eventualbuddha/zedbar/compare/v0.5.0...v0.5.1) (2026-08-15)


### Bug Fixes

* **cli:** say why an image could not be scanned ([792961a](https://github.com/eventualbuddha/zedbar/commit/792961a1638ad8218889041da1d597c915e1a576))
* **qrcode:** correct GS1, encoding and structured-append decoding, and five overflow traps ([daa7edd](https://github.com/eventualbuddha/zedbar/commit/daa7eddd974285122142bcb3c1ab2648c6e2536a))
* **qrcode:** correct three faults in QR text extraction ([ce40222](https://github.com/eventualbuddha/zedbar/commit/ce40222a045425e70e5bc7ff47a1872b2c14cb07))
* **qrcode:** keep the alignment-pattern centre inside the image ([932002c](https://github.com/eventualbuddha/zedbar/commit/932002c8c2e55f3a7fa8ab5f9f6b7c10f707b24b))
* **qrcode:** make the second homography builder survive runaway corners ([34ab80e](https://github.com/eventualbuddha/zedbar/commit/34ab80e496977be05c45c7d73ca0fd548678b5b4))
* **qrcode:** reject an alignment-pattern corner outside the image ([08e223e](https://github.com/eventualbuddha/zedbar/commit/08e223e99c2d3ccd49af927d6801a43e99836ace))
* **qrcode:** report an incomplete structured-append group as partial ([d503fdc](https://github.com/eventualbuddha/zedbar/commit/d503fdccc77b7cc82f619f2bd2695f7f71128550))
* **qrcode:** saturate the turn test the way the distance test already does ([c8bcea5](https://github.com/eventualbuddha/zedbar/commit/c8bcea54cb8656394329eb812c1f6c329a206eef))
* **qrcode:** stop three more corner computations from trapping ([c29b7a5](https://github.com/eventualbuddha/zedbar/commit/c29b7a570903c7e3f1afab408c88841353425695))
* reject a zero-width reference character instead of dividing by it ([cc927f2](https://github.com/eventualbuddha/zedbar/commit/cc927f29317af5b8c8ae5157c1cffa82ed7a8e44))
* **scanner:** take only 2D symbols from a retried finder region ([14f9ce0](https://github.com/eventualbuddha/zedbar/commit/14f9ce03af27298d899ef69b5391748a4272ad13))
* **zedbarimg:** keep incomplete symbols out of the decoded output ([018449b](https://github.com/eventualbuddha/zedbar/commit/018449bc5c7182da4a036053aa0289ced17f2e43))

## [0.5.0](https://github.com/eventualbuddha/zedbar/compare/v0.4.2...v0.5.0) (2026-08-15)


### ⚠ BREAKING CHANGES

* `Symbol::data()` returns the raw bit matrix for SQ Code symbols rather than its base64 encoding, and `Symbol::raw_data()` returns `None` for them. Callers wanting the previous value should base64-encode `data()`. `zedbarimg` prints the decoded bytes where it used to print base64.

### Bug Fixes

* correct EAN add-on indexing, checksum-suppressed length, and i25 exit check ([e174c8c](https://github.com/eventualbuddha/zedbar/commit/e174c8c1613ac2ff726e874c5a22d8fc608d7ee0))
* **databar:** repair the output buffer layout and expanded bit feeder ([e32104f](https://github.com/eventualbuddha/zedbar/commit/e32104faa37931c46da6b89118bdfc7348cbd5b0))
* **qrcode:** bound the undecoded-region search and the retry that consumes it ([4f71022](https://github.com/eventualbuddha/zedbar/commit/4f710227f4269474e534876cc50b9c9d000acb84))
* **qrcode:** read the corner bounds shift with C's operator precedence ([49663a9](https://github.com/eventualbuddha/zedbar/commit/49663a9d59b6033e7064bab132cb95ec5bfa6d6d))
* report SQ Code payloads as bytes instead of base64 ([0134c37](https://github.com/eventualbuddha/zedbar/commit/0134c37a7de0a956c82f732956aa34610d596911))
* reset decoder state between images and cap regions on every exit ([9b8a49f](https://github.com/eventualbuddha/zedbar/commit/9b8a49fa30ec25dafdb7515c9cc8156b6b34cdb6))
* **scanner:** reject flat 2nd differentials when detecting edges ([044a8d5](https://github.com/eventualbuddha/zedbar/commit/044a8d56a31e4ecefcba80237d322cd5975eb948))
* **sqcode:** stop reporting the base64 payload with a trailing NUL ([74277c4](https://github.com/eventualbuddha/zedbar/commit/74277c426c929db93d529262ed0821a18017c09b))
* stop imposing a four-character minimum on Code 93 and Code 128 ([0659320](https://github.com/eventualbuddha/zedbar/commit/06593204bbe3b24c53cc1ed6eca0b95736a557c6))

## [0.4.2](https://github.com/eventualbuddha/zedbar/compare/v0.4.1...v0.4.2) (2026-08-13)


### Bug Fixes

* **codabar:** compute width products in 64 bits like the C decoder ([d58f3e0](https://github.com/eventualbuddha/zedbar/commit/d58f3e039bce5392bb2bfd8de52c31474f60dbf9))
* **code128:** preserve character count during postprocess ([#61](https://github.com/eventualbuddha/zedbar/issues/61)) ([02f9258](https://github.com/eventualbuddha/zedbar/commit/02f92580942435b8a351f3b8c69985e865a9c919))
* **code39:** report the stop character as '?' and bound the buffer ([a759d67](https://github.com/eventualbuddha/zedbar/commit/a759d67536c4e52d1f6d5dae71bbc1eb55a5ac09))
* **code93:** take the Code 93 lock and abort on a rejected stop ([b791107](https://github.com/eventualbuddha/zedbar/commit/b791107ba1ba2c44982b487aa84b61344fcd9a48))
* **decoder:** only release the shared lock to the symbology holding it ([e035b08](https://github.com/eventualbuddha/zedbar/commit/e035b089b7538a2a7def65208cc77bc70399e8e3))
* **demo:** fill the canvas white before drawing the image ([7eebd27](https://github.com/eventualbuddha/zedbar/commit/7eebd271dfe41465a6ca2f2188c0b104b590e47c))
* **ean:** stop panicking when merging incompatible partial halves ([c3a0e73](https://github.com/eventualbuddha/zedbar/commit/c3a0e73704fbb275b9a326a2b1be00eb09eef2ec))
* **i25:** reject invalid element widths and bound the holding buffer ([b38e9a5](https://github.com/eventualbuddha/zedbar/commit/b38e9a5b26402d0c6f11766a1cf006aa9d9d1a9f))
* **image:** composite transparency over white before scanning ([6931f43](https://github.com/eventualbuddha/zedbar/commit/6931f43dcb04423bf9d4b44398cf225f1b9e721c))
* **image:** stop panicking on empty images and out-of-range crops ([3c7bde4](https://github.com/eventualbuddha/zedbar/commit/3c7bde4c6a42b830a24683bfb4b8cf9b3d947459))
* **symbol:** report Right and Down the way zbar defines them ([199989b](https://github.com/eventualbuddha/zedbar/commit/199989b5b723b0a49b276ed9dc8c2fc5f2013c5e))

## [0.4.1](https://github.com/eventualbuddha/zedbar/compare/v0.4.0...v0.4.1) (2026-05-11)


### Bug Fixes

* **qrcode:** report tight per-triplet failure regions for retry ([#55](https://github.com/eventualbuddha/zedbar/issues/55)) ([f80dcb6](https://github.com/eventualbuddha/zedbar/commit/f80dcb685aad8c9972d4515ac0669ec1fba959b1))

## [0.4.0](https://github.com/eventualbuddha/zedbar/compare/v0.3.1...v0.4.0) (2026-05-08)


### Features

* expose symbol position points and bounding rectangle ([#52](https://github.com/eventualbuddha/zedbar/issues/52)) ([4ab6720](https://github.com/eventualbuddha/zedbar/commit/4ab67200aadb5a9ed0e4bc599b975ee814b94367))

## [0.3.1](https://github.com/eventualbuddha/zedbar/compare/v0.3.0...v0.3.1) (2026-05-07)


### Bug Fixes

* better handle noisy images ([#50](https://github.com/eventualbuddha/zedbar/issues/50)) ([1d5db4a](https://github.com/eventualbuddha/zedbar/commit/1d5db4a3e53515da4633a661785bbe554eabdd8d))

## [0.3.0](https://github.com/eventualbuddha/zedbar/compare/v0.2.5...v0.3.0) (2026-05-06)


### Features

* make `DecoderConfig::new()` start empty ([#48](https://github.com/eventualbuddha/zedbar/issues/48)) ([269a961](https://github.com/eventualbuddha/zedbar/commit/269a9614f2e356cb055692296ff9c2ec31edab2c))

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
