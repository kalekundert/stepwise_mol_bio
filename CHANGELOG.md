# Changelog

<!--next-version-placeholder-->

## v1.23.1 (2021-06-29)
### Fix
* **laser_scanner:** Correctly use list-based presets ([`c8fc44f`](https://github.com/kalekundert/stepwise_mol_bio/commit/c8fc44fcbc1eea9dde196444dc9d730ed306ce91))

## v1.23.0 (2021-06-28)
### Feature
* **lyophilize:** Add protocol ([`4c43fb0`](https://github.com/kalekundert/stepwise_mol_bio/commit/4c43fb01e5e4962ef0ae5e02af4cc688471ed0fe))
* Upgrade appcli; get protocols from products via CLI ([`fa19df2`](https://github.com/kalekundert/stepwise_mol_bio/commit/fa19df29e1ca207b29a208b508c63f85385d9ec8))
* **digest:** Report product concentration/volume ([`b04ae06`](https://github.com/kalekundert/stepwise_mol_bio/commit/b04ae06928b9a24f03e9f2052d775ec9313f5faa))
* **miniprep:** Add very simple protocol ([`7a15c32`](https://github.com/kalekundert/stepwise_mol_bio/commit/7a15c32053107be139e0dc2a6515a6b9105f5c95))

### Documentation
* Fix typo ([`020bd3e`](https://github.com/kalekundert/stepwise_mol_bio/commit/020bd3ec5ddabcd4bd3613e12e3e6fd8651b5136))

## v1.22.0 (2021-06-08)
### Feature
* Tweak wording ([`f53013d`](https://github.com/kalekundert/stepwise_mol_bio/commit/f53013dc3c3bdbabd4f793c68989aec6fb9426d8))
* **ligate:** Implement freezerbox maker plugin ([`5a98321`](https://github.com/kalekundert/stepwise_mol_bio/commit/5a98321cb57b23decb76de72e9bd654ce9b7ae50))
* **gibson:** Implement freezerbox maker plugin ([`ec5df7d`](https://github.com/kalekundert/stepwise_mol_bio/commit/ec5df7da03d3fc148b1df0f5663d4b7f57ad372e))

## v1.21.0 (2021-06-07)
### Feature
* **golden_gate:** Implement freezerbox maker plugin ([`cd53b3e`](https://github.com/kalekundert/stepwise_mol_bio/commit/cd53b3e1f28b4bad105f79e8ce0d0c9afa05ec0b))

## v1.20.0 (2021-06-02)
### Feature
* **page_purify:** Add protocol ([`ba10309`](https://github.com/kalekundert/stepwise_mol_bio/commit/ba103099dea789390e57b3eea1880ae6d8783578))
* **gel:** Add pre-running options and config docs ([`4911a75`](https://github.com/kalekundert/stepwise_mol_bio/commit/4911a75f548685337f7f9d47b647fd53c586e517))

## v1.19.0 (2021-05-19)
### Feature
* **spin_cleanup:** Add Zymo RNA Clean & Concentrator protocols ([`713f444`](https://github.com/kalekundert/stepwise_mol_bio/commit/713f444164c4de97ea8605cd94e64869b1139e1d))

### Fix
* **ivt:** Update variable name ([`c3d73ce`](https://github.com/kalekundert/stepwise_mol_bio/commit/c3d73ce43f057888e3715943d5304ce683cfa328))

## v1.18.0 (2021-05-11)
### Feature
* **spin_cleanup:** Implement freezerbox maker plugin ([`b77cda9`](https://github.com/kalekundert/stepwise_mol_bio/commit/b77cda92fa8991c30b21dbe593e3987459f2029e))

### Fix
* **gel:** Don't try to adjust volumes if there is no solvent ([`cb7f4af`](https://github.com/kalekundert/stepwise_mol_bio/commit/cb7f4af8f493bb800d77051d7260b3ca38714adb))

## v1.17.0 (2021-05-05)
### Feature
* **digest:** Implement freezerbox maker plugin ([`b95a73f`](https://github.com/kalekundert/stepwise_mol_bio/commit/b95a73ffde8f94c8c6312dc5eb199351e4a944b7))
* **ivtt:** Add an --extra-percent option ([`5e25d56`](https://github.com/kalekundert/stepwise_mol_bio/commit/5e25d56519019f5228e885df8050128aca9e32de))
* **aliquot:** Implement freezerbox maker plugin ([`4b7977f`](https://github.com/kalekundert/stepwise_mol_bio/commit/4b7977fe6445c0613de5ed8cffac3f85efbbdefb))

### Fix
* **ivt:** Include rNTP mix in the master mix ([`6306ee6`](https://github.com/kalekundert/stepwise_mol_bio/commit/6306ee628acaa5ba70fb13b4ee1d1376eec107eb))
* **spin_cleanup:** Update preset name ([`2eb8419`](https://github.com/kalekundert/stepwise_mol_bio/commit/2eb8419a8a816477ffe0ab06fa8db8ca8bba9f7c))
* **ivt:** Don't require unanimous template concentrations ([`1496226`](https://github.com/kalekundert/stepwise_mol_bio/commit/14962261341bd414cb5994e9f4dbb93f9d8e41e8))

### Documentation
* **ivt:** Tweak database section ([`9bbb70f`](https://github.com/kalekundert/stepwise_mol_bio/commit/9bbb70f219e21e145fdb393c617e55cdaaf63e23))

## v1.16.0 (2021-04-26)
### Feature
* **ivt:** Implement freezerbox maker plugin ([`b07cdfb`](https://github.com/kalekundert/stepwise_mol_bio/commit/b07cdfbc65f97d5edc9e513c8cf7f308bdeb5aa7))
* **rnasezap:** Add protocol ([`6ea3401`](https://github.com/kalekundert/stepwise_mol_bio/commit/6ea3401f39a4c07a97b54fb4b72811f23d1bb6f8))
* **invpcr:** Implement freezerbox maker plugin ([`fab5240`](https://github.com/kalekundert/stepwise_mol_bio/commit/fab5240ccc28c44fb93681d5297a2cc8dc36fc2a))

### Documentation
* Remove docs badge ([`c2bd23d`](https://github.com/kalekundert/stepwise_mol_bio/commit/c2bd23ddd6c99a4e564aa6250b4a8007060cd004))

## v1.15.0 (2021-04-23)
### Feature
* **pcr:** Implement freezerbox maker plugin ([`5cf472a`](https://github.com/kalekundert/stepwise_mol_bio/commit/5cf472a1178ea5f61c6127654906f9f84e1e2e22))
* **gel:** Add tricine gel parameters ([`240c749`](https://github.com/kalekundert/stepwise_mol_bio/commit/240c7493dfcb42c5bc53a5d6fb1e2764d367e91b))
* **spin_cleanup:** Add a sample dilution step ([`99868f6`](https://github.com/kalekundert/stepwise_mol_bio/commit/99868f6d62dfd1ccc48b97349216058b587115d3))
* **ivt:** Add --volume-uL option ([`109b0a9`](https://github.com/kalekundert/stepwise_mol_bio/commit/109b0a95d4ac31869b8619cdf485bf573e2d5c22))
* Add generic staining protocol presets ([`4df2876`](https://github.com/kalekundert/stepwise_mol_bio/commit/4df2876c1f8e2c6ca67c018b3fb81227f734bb39))

### Fix
* **stain:** Debug -I flag ([`9cda251`](https://github.com/kalekundert/stepwise_mol_bio/commit/9cda2519193b0adc68e0e4e579d807ba446c7f43))

## v1.14.0 (2021-04-02)
### Feature
* **ivtt:** Use exact times, not ranges ([`0f01ef1`](https://github.com/kalekundert/stepwise_mol_bio/commit/0f01ef1383906aea13a2bc3606c00063427eecba))
* **gel:** Make preset names more consistent ([`a824515`](https://github.com/kalekundert/stepwise_mol_bio/commit/a82451584b93aebc8d0904c0723e5f6ed3d02c0c))

### Fix
* **ivtt:** Strip whitespace when parsing reagents ([`c509e03`](https://github.com/kalekundert/stepwise_mol_bio/commit/c509e0374f5f771d55b93c80591a27a33bb6cf26))
* **gel:** Tweak wording ([`d716e73`](https://github.com/kalekundert/stepwise_mol_bio/commit/d716e733eae2fd8fcf59d2071d7bdae192b34196))
* **ivt:** Correct T7 volume for short templates ([`f41f74b`](https://github.com/kalekundert/stepwise_mol_bio/commit/f41f74b4b2a4f45edaaf63a8d18b35ff712dbb5e))

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
