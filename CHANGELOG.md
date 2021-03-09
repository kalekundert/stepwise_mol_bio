# Changelog

<!--next-version-placeholder-->

## v1.13.0 (2021-03-09)
### Feature
* **ivtt:** Add a PUREfrex 1.0 protocol ([`8d87531`](https://github.com/kalekundert/stepwise_mol_bio/commit/8d87531209e6417bb29784bb1bfb2a0691882d49))

### Fix
* **serial:** Update formatting ([`e01791e`](https://github.com/kalekundert/stepwise_mol_bio/commit/e01791e035cd39ebe14c73c55b6ef59abf12d258))

## v1.12.0 (2021-03-09)
### Feature
* **spin_column:** Include title in protocol ([`58235a3`](https://github.com/kalekundert/stepwise_mol_bio/commit/58235a3d36791cb35904157134335b5bf1d47da6))

### Fix
* **pcr:** Update formatting ([`4747f32`](https://github.com/kalekundert/stepwise_mol_bio/commit/4747f328833d66d2ad7b337d8c840d7a88c42c9c))
* Remove imports of deleted modules ([`e170378`](https://github.com/kalekundert/stepwise_mol_bio/commit/e170378600ac1a654adb1c2902929b1633669a19))
* **gel:** Update formatting ([`63b65e2`](https://github.com/kalekundert/stepwise_mol_bio/commit/63b65e209ef6d89980fce71672aad0e62ce11007))

## v1.11.0 (2021-03-08)
### Feature
* **spin_cleanup:** Add protocol ([`1eda6f9`](https://github.com/kalekundert/stepwise_mol_bio/commit/1eda6f98074a0c659be74e1fb3e974439cfccf2d))
* **stain:** Add a general staining protocol ([`83767be`](https://github.com/kalekundert/stepwise_mol_bio/commit/83767be74410bbbc0614cfe179998a88b09c5d89))

### Fix
* Use `ul()`/`pl()` instead of `Step()` ([`af42c54`](https://github.com/kalekundert/stepwise_mol_bio/commit/af42c5440bfe4d2a36fd927dccdc8df4813c0402))

## v1.10.0 (2021-03-01)
### Feature
* Be more explicit about when to add nuclease inhibitors ([`eee734e`](https://github.com/kalekundert/stepwise_mol_bio/commit/eee734e68538808c99fe821af53c2a9532a46eaa))

### Fix
* Interpret additive volumes in the right context ([`855108a`](https://github.com/kalekundert/stepwise_mol_bio/commit/855108aff13032402b05e4009386d8e3938648bd))

## v1.9.1 (2021-02-15)
### Fix
* **ethanol:** Use correct buffer volume option ([`4135336`](https://github.com/kalekundert/stepwise_mol_bio/commit/4135336d583368e84ece8c16fc385451ccc3ba1b))

## v1.9.0 (2021-02-14)
### Feature
* Add the `ivtt` protocol ([`b3426d4`](https://github.com/kalekundert/stepwise_mol_bio/commit/b3426d49e8e52ebd864f5c861d5ae93ce9a5cea4))

## v1.8.0 (2021-02-03)
### Feature
* **gel:** Include buffer in protocol ([`3fabfb8`](https://github.com/kalekundert/stepwise_mol_bio/commit/3fabfb82fc91febac59687d9c6ac796f0d012fb7))
* Add option to skip cleanup step ([`95e4724`](https://github.com/kalekundert/stepwise_mol_bio/commit/95e472403e17dbbe1eca3a8f6be8cc6bb2867878))

## v1.7.0 (2021-01-20)
### Feature
* Add option to override template stock conc ([`8280110`](https://github.com/kalekundert/stepwise_mol_bio/commit/8280110228d69cb5a0b496aee4c9a20ed78962d1))

## v1.6.0 (2021-01-11)
### Feature
* Migrate to appcli ([`342f863`](https://github.com/kalekundert/stepwise_mol_bio/commit/342f8637cac8b5ac1d36ac0d9f1f19c6db883cc6))

## v1.5.0 (2020-11-03)
### Feature
* Use fresh 70% ethanol ([`eb44a92`](https://github.com/kalekundert/stepwise_mol_bio/commit/eb44a92cc6f947a87343cee03cf8116e531e7897))
* Implement ethanol precipitation protocol from Li2020 ([`a6e79ff`](https://github.com/kalekundert/stepwise_mol_bio/commit/a6e79ffb2ea9683a8df40f3558f71e0363caaa1b))

### Fix
* Reorganize footnotes for ethanol precipitation ([`68a24b9`](https://github.com/kalekundert/stepwise_mol_bio/commit/68a24b9c046847be04657e643e859297932c21be))

## v1.4.1 (2020-10-20)
### Fix
* Put the phenol-chloroform protocol in the right directory ([`934c07d`](https://github.com/kalekundert/stepwise_mol_bio/commit/934c07daaf56d5ea8a96100e1b3d08d84b6ddca0))

## v1.4.0 (2020-10-20)
### Feature
* Add a phenol-chloroform extraction protocol ([`fe632e1`](https://github.com/kalekundert/stepwise_mol_bio/commit/fe632e1c1bcfac0ab33841f9cc381e6ead0556e9))

## v1.3.0 (2020-10-19)
### Feature
* Add the aliquot protocol ([`b79b206`](https://github.com/kalekundert/stepwise_mol_bio/commit/b79b2066eee246beeb4fd4868623b74445465ee1))
* Add flags to skip either step of inverse PCR ([`07bdf44`](https://github.com/kalekundert/stepwise_mol_bio/commit/07bdf44a6023e2683eb4786be8bfec36e19b69fe))
* Add the -n flag for restriction digests ([`e80f0f8`](https://github.com/kalekundert/stepwise_mol_bio/commit/e80f0f80f44994b95cc0d558b35a26b51c6186e8))
* Add flag for gradient PCR ([`f7e9454`](https://github.com/kalekundert/stepwise_mol_bio/commit/f7e945469b38a355bb01253f38c8c51d2bd64ff8))
* Make imaging step optional (for staining protocols) ([`2cf793f`](https://github.com/kalekundert/stepwise_mol_bio/commit/2cf793fdff6359fee45aead23a5458b77683a517))

### Fix
* Make the time parameters integers ([`165dca9`](https://github.com/kalekundert/stepwise_mol_bio/commit/165dca9cabf46dfb3d0ed14eec3ee53e6cda6622))
* Give salt concentrations in molarity instead of percent ([`15a4693`](https://github.com/kalekundert/stepwise_mol_bio/commit/15a46938cf01d6f8851bbc72bbd69f6474a5e1a9))
* Correct primer volumes for scaled PCR reactions ([`3bc8d42`](https://github.com/kalekundert/stepwise_mol_bio/commit/3bc8d422ea8bfb9a5cfcbef79d170174f90b34be))
