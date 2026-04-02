# Changelog

## [0.4.1](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.4.0...v0.4.1) (2026-04-02)


### Features

* export in DNPM 2.1 JSON format ([#30](https://github.com/pcvolkmer/mhguide-to-csv/issues/30)) ([ce3f95e](https://github.com/pcvolkmer/mhguide-to-csv/commit/ce3f95eedf2fea5148dfada738207df00044a32e))
* export read depth for simple variants ([#31](https://github.com/pcvolkmer/mhguide-to-csv/issues/31)) ([aab6022](https://github.com/pcvolkmer/mhguide-to-csv/commit/aab6022e86c9e4ba80c810c283e79a081b80cb83))

## [0.4.0](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.9...v0.4.0) (2026-03-26)


### ⚠ BREAKING CHANGES

* split records in multiple parts ([#28](https://github.com/pcvolkmer/mhguide-to-csv/issues/28))

### Features

* convert msi value to % ([8b2bfde](https://github.com/pcvolkmer/mhguide-to-csv/commit/8b2bfde6586808f45d007f7afb190d8a86e9d910))
* export RNA fusions ([#29](https://github.com/pcvolkmer/mhguide-to-csv/issues/29)) ([bfc2d7a](https://github.com/pcvolkmer/mhguide-to-csv/commit/bfc2d7acfc3db4feb806610b7944d435bd097d3e))
* parse RNA fusions from string input ([#26](https://github.com/pcvolkmer/mhguide-to-csv/issues/26)) ([8720f9f](https://github.com/pcvolkmer/mhguide-to-csv/commit/8720f9f15d16023e1388dab3a8f369df93c2444f))
* split records in multiple parts ([#28](https://github.com/pcvolkmer/mhguide-to-csv/issues/28)) ([51727f2](https://github.com/pcvolkmer/mhguide-to-csv/commit/51727f2f0b70eabfec6b8105dbadecc8bdd3b2ac))


### Bug Fixes

* biomarkers format to use comma-separated decimals ([00a3700](https://github.com/pcvolkmer/mhguide-to-csv/commit/00a3700bc22a101d83a1cda833ae6f918c3f1fb6))
* use msi score value as is ([837bdae](https://github.com/pcvolkmer/mhguide-to-csv/commit/837bdae2e8d345bb491eff1374803e69ea422981))

## [0.3.9](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.8...v0.3.9) (2026-03-16)


### Features

* add DISPLAY_MODIFIED_OBJECT to biomarkers ([c948439](https://github.com/pcvolkmer/mhguide-to-csv/commit/c9484394c3d92e64214b1b3f7e518b792325f987))
* add parsing of biomarkers TMB, HRD and MSI ([2b52ef5](https://github.com/pcvolkmer/mhguide-to-csv/commit/2b52ef5ca1eb01419cb857eb06c16cf3286e9ded))
* export HRD, MSI and TMB ([#25](https://github.com/pcvolkmer/mhguide-to-csv/issues/25)) ([ebd7445](https://github.com/pcvolkmer/mhguide-to-csv/commit/ebd74459669506101f4b3b1981128d17d9029e54))


### Bug Fixes

* use TMB variant type ([e93d5ac](https://github.com/pcvolkmer/mhguide-to-csv/commit/e93d5ac5a108bfcf0bc9422fe97473f99a9e7807))

## [0.3.8](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.7...v0.3.8) (2026-03-07)


### Features

* make artifact removing optional ([8630074](https://github.com/pcvolkmer/mhguide-to-csv/commit/8630074267e1c9c960500c8e83495a6ae2ebe33a))

## [0.3.7](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.6...v0.3.7) (2026-03-06)


### Features

* add biomarker structs ([38afacb](https://github.com/pcvolkmer/mhguide-to-csv/commit/38afacb983e07479758a9d97fd80d5c586fd2ed9))
* add cnv from biomarkers ([51ae04e](https://github.com/pcvolkmer/mhguide-to-csv/commit/51ae04eccd800d29eec00b7a9e850f44ffcc8bbe))
* add more CNV from report narrative ([128db30](https://github.com/pcvolkmer/mhguide-to-csv/commit/128db30fd7da816dffe8a5314441d4c1a77e72cf))

## [0.3.6](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.5...v0.3.6) (2026-03-04)


### Bug Fixes

* add other used artifact wordings ([42d77fb](https://github.com/pcvolkmer/mhguide-to-csv/commit/42d77fbbac00dd168e8381b3469c4740ca749a63))

## [0.3.5](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.4...v0.3.5) (2026-03-02)


### Features

* add CNV from report narrative ([#18](https://github.com/pcvolkmer/mhguide-to-csv/issues/18)) ([abeb67a](https://github.com/pcvolkmer/mhguide-to-csv/commit/abeb67a4055324fbffb138fd9b5dd7340f702dd1))

## [0.3.4](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.3...v0.3.4) (2026-03-02)


### Features

* check zip file content for single JSON file ([e35369c](https://github.com/pcvolkmer/mhguide-to-csv/commit/e35369c8e235ad5ceb67627e4d662ca9720abff3))
* remove artifact variants ([a4d06c8](https://github.com/pcvolkmer/mhguide-to-csv/commit/a4d06c8d9ac53463ccc61d55e0fa784219ee648a))

## [0.3.3](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.2...v0.3.3) (2026-02-26)


### Features

* accept ZIP files ([#16](https://github.com/pcvolkmer/mhguide-to-csv/issues/16)) ([7f6fb68](https://github.com/pcvolkmer/mhguide-to-csv/commit/7f6fb68225469eee845a17b448f98eb1f2292232))
* export into Excel file ([#14](https://github.com/pcvolkmer/mhguide-to-csv/issues/14)) ([00b68e9](https://github.com/pcvolkmer/mhguide-to-csv/commit/00b68e91777e3a7451b91b70bb804eb0086c5865))

## [0.3.2](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.1...v0.3.2) (2026-02-26)


### Features

* add pathogenic classification ([#13](https://github.com/pcvolkmer/mhguide-to-csv/issues/13)) ([ec45c19](https://github.com/pcvolkmer/mhguide-to-csv/commit/ec45c197435aa60d64a00a01612a6d18ec5413d9))
* search genes by previous symbols ([#12](https://github.com/pcvolkmer/mhguide-to-csv/issues/12)) ([61e27b8](https://github.com/pcvolkmer/mhguide-to-csv/commit/61e27b88d3d3620d178cec8b188a4c6dde1d5cc3))


### Bug Fixes

* allow hyphen in gene symbol ([dadd3fa](https://github.com/pcvolkmer/mhguide-to-csv/commit/dadd3faa052ae64d8b5a593ab01275ea1131d396))
* allow underscore in gene symbol ([7f5406d](https://github.com/pcvolkmer/mhguide-to-csv/commit/7f5406d1eed046ccc71cb545a9e52ac266d54ed1))

## [0.3.1](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.3.0...v0.3.1) (2026-02-25)


### Features

* add narrative report variants with cDNA modification ([#9](https://github.com/pcvolkmer/mhguide-to-csv/issues/9)) ([7adf573](https://github.com/pcvolkmer/mhguide-to-csv/commit/7adf5731b1d169ef4ba1f30f48bafa33b6191983))

## [0.3.0](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.2.2...v0.3.0) (2026-02-25)


### ⚠ BREAKING CHANGES

* include variants from `REPORT_NARRATIVE` ([#7](https://github.com/pcvolkmer/mhguide-to-csv/issues/7))

### Features

* include variants from `REPORT_NARRATIVE` ([#7](https://github.com/pcvolkmer/mhguide-to-csv/issues/7)) ([f0409c1](https://github.com/pcvolkmer/mhguide-to-csv/commit/f0409c1e3aaf669f635912cedf64dad99068a0e7))

## [0.2.2](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.2.1...v0.2.2) (2026-02-19)


### Features

* add arg to use oncogenic classification variants only ([a4dd8df](https://github.com/pcvolkmer/mhguide-to-csv/commit/a4dd8df2b453ba10edf24f9c3099b1d019df834b))

## [0.2.1](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.2.0...v0.2.1) (2026-02-17)


### Bug Fixes

* include all benign oncogenic classification variants ([#4](https://github.com/pcvolkmer/mhguide-to-csv/issues/4)) ([88b0199](https://github.com/pcvolkmer/mhguide-to-csv/commit/88b0199b7eba89938d8529fc62035a8b2fee7f5b))

## [0.2.0](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.1.0...v0.2.0) (2026-02-17)


### ⚠ BREAKING CHANGES

* filter required variants based on oncogenic classification ([#2](https://github.com/pcvolkmer/mhguide-to-csv/issues/2))

### Features

* filter required variants based on oncogenic classification ([#2](https://github.com/pcvolkmer/mhguide-to-csv/issues/2)) ([2f0866e](https://github.com/pcvolkmer/mhguide-to-csv/commit/2f0866ee056e00be12508892ec87ba8ad1d8cc50))

## [0.1.0](https://github.com/pcvolkmer/mhguide-to-csv/compare/v0.1.0...v0.1.0) (2026-02-17)


### Features

* add duplications ([aefa21d](https://github.com/pcvolkmer/mhguide-to-csv/commit/aefa21d1dd521c5630d1fbd04a9fe055af598c16))
* add hgnc gene list to complete result list ([2e2a859](https://github.com/pcvolkmer/mhguide-to-csv/commit/2e2a8590f3c3821512d42b6f780a9c416b7e9d40))
* map protein modifications to three-letter format ([06dd3ef](https://github.com/pcvolkmer/mhguide-to-csv/commit/06dd3ef9670342e62a54147d9ab2451716d6006d))
* show original protein modification ([7277db7](https://github.com/pcvolkmer/mhguide-to-csv/commit/7277db73636ab2c5d778cf81b51eeb8e61a64f0a))
* support multi insertions ([c1fc8ae](https://github.com/pcvolkmer/mhguide-to-csv/commit/c1fc8ae206e856f942a92a406f74fb10c69e86b7))


### Bug Fixes

* p. insdel three letter code ([662d3df](https://github.com/pcvolkmer/mhguide-to-csv/commit/662d3dffe60809643a15b269528e1b143c662f92))
* use chromosome_modification for start/end ([3dcaa9e](https://github.com/pcvolkmer/mhguide-to-csv/commit/3dcaa9eace774ed17bac9538f8b73eeb21179fda))


### Miscellaneous Chores

* release 0.1.0 ([539c7bb](https://github.com/pcvolkmer/mhguide-to-csv/commit/539c7bb8840d70c55eee57a23ff19c5ede5aa4f6))
