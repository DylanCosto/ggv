# Third-party notices

GGV is MIT-licensed (see [LICENSE](LICENSE)). Its release binaries (`ggv`, `ggv-desktop`) include the
open-source crates below (191 in all), each under its own license as its authors publish it.
Every license here is permissive; where a crate offers a choice (for example `MIT OR Apache-2.0`, or
`Apache-2.0 OR GPL-2.0-only`), GGV uses it under the permissive option. The full license texts are in
each crate's repository and in its package on crates.io.

Fonts: the interface uses the fonts bundled with egui (`epaint_default_fonts`): Ubuntu (Ubuntu Font
Licence 1.0), Hack (MIT, with Bitstream Vera and DejaVu notices), Noto Emoji (SIL Open Font License
1.1) and emoji-icon-font (MIT or Apache-2.0).

Data: GGV embeds a compact copy of the MANE GRCh38 v1.5 summary (Matched Annotation from NCBI and
EMBL-EBI; RefSeq and Ensembl transcript ids with their MANE Select / MANE Plus Clinical status, in
`crates/ggv-annot/data/mane-grch38-v1.5.tsv.gz`), used to pick each gene's main transcript. MANE is
a public-domain resource from NCBI and EMBL-EBI (https://www.ncbi.nlm.nih.gov/refseq/MANE/).

Regenerate this file with `cargo tree -e normal -p ggv-desktop -p ggv-cli --prefix none --format "{p}|{l}|{r}"`
whenever dependencies change.

## Licenses at a glance

| License | Crates |
|---|---|
| MIT OR Apache-2.0 | 113 |
| Apache-2.0 OR MIT | 18 |
| MIT | 14 |
| MIT/Apache-2.0 | 9 |
| Apache-2.0 | 7 |
| MIT OR Apache-2.0 OR Zlib | 3 |
| (MIT OR Apache-2.0) AND Unicode-3.0 | 2 |
| MIT OR Zlib OR Apache-2.0 | 2 |
| Zlib | 2 |
| Unlicense OR MIT | 2 |
| BSD-3-Clause | 2 |
| Zlib OR Apache-2.0 OR MIT | 2 |
| BSD-2-Clause OR Apache-2.0 OR MIT | 2 |
| Apache-2.0/MIT | 2 |
| BSL-1.0 | 2 |
| BSD-3-Clause OR Apache-2.0 | 2 |
| 0BSD OR MIT OR Apache-2.0 | 1 |
| bzip2-1.0.6 | 1 |
| BSD-2-Clause | 1 |
| (MIT OR Apache-2.0) AND OFL-1.1 AND Ubuntu-font-1.0 | 1 |
| Apache-2.0 OR GPL-2.0-only | 1 |
| ISC | 1 |
| Apache-2.0 AND MIT | 1 |

## Crates

| Crate | Version | License | Source |
|---|---|---|---|
| accesskit | 0.25.2 | MIT OR Apache-2.0 | [github.com/AccessKit/accesskit](https://github.com/AccessKit/accesskit) |
| accesskit_consumer | 0.39.2 | MIT OR Apache-2.0 | [github.com/AccessKit/accesskit](https://github.com/AccessKit/accesskit) |
| accesskit_windows | 0.35.2 | MIT OR Apache-2.0 | [github.com/AccessKit/accesskit](https://github.com/AccessKit/accesskit) |
| accesskit_winit | 0.34.2 | Apache-2.0 | [github.com/AccessKit/accesskit](https://github.com/AccessKit/accesskit) |
| adler2 | 2.0.1 | 0BSD OR MIT OR Apache-2.0 | [github.com/oyvindln/adler2](https://github.com/oyvindln/adler2) |
| ahash | 0.8.12 | MIT OR Apache-2.0 | [github.com/tkaitchuck/ahash](https://github.com/tkaitchuck/ahash) |
| allocator-api2 | 0.2.21 | MIT OR Apache-2.0 | [github.com/zakarumych/allocator-api2](https://github.com/zakarumych/allocator-api2) |
| anstream | 1.0.0 | MIT OR Apache-2.0 | [github.com/rust-cli/anstyle.git](https://github.com/rust-cli/anstyle.git) |
| anstyle | 1.0.14 | MIT OR Apache-2.0 | [github.com/rust-cli/anstyle.git](https://github.com/rust-cli/anstyle.git) |
| anstyle-parse | 1.0.0 | MIT OR Apache-2.0 | [github.com/rust-cli/anstyle.git](https://github.com/rust-cli/anstyle.git) |
| anstyle-query | 1.1.5 | MIT OR Apache-2.0 | [github.com/rust-cli/anstyle.git](https://github.com/rust-cli/anstyle.git) |
| anstyle-wincon | 3.0.11 | MIT OR Apache-2.0 | [github.com/rust-cli/anstyle.git](https://github.com/rust-cli/anstyle.git) |
| anyhow | 1.0.104 | MIT OR Apache-2.0 | [github.com/dtolnay/anyhow](https://github.com/dtolnay/anyhow) |
| arboard | 3.7.0 | MIT OR Apache-2.0 | [github.com/1Password/arboard](https://github.com/1Password/arboard) |
| arrayref | 0.3.9 | BSD-2-Clause | [github.com/droundy/arrayref](https://github.com/droundy/arrayref) |
| arrayvec | 0.7.8 | MIT OR Apache-2.0 | [github.com/bluss/arrayvec](https://github.com/bluss/arrayvec) |
| ash | 0.38.0+1.3.281 | MIT OR Apache-2.0 | [github.com/ash-rs/ash](https://github.com/ash-rs/ash) |
| bit-set | 0.10.0 | Apache-2.0 OR MIT | [github.com/contain-rs/bit-set](https://github.com/contain-rs/bit-set) |
| bit-vec | 0.10.1 | Apache-2.0 OR MIT | [github.com/contain-rs/bit-vec](https://github.com/contain-rs/bit-vec) |
| bit-vec | 0.9.1 | Apache-2.0 OR MIT | [github.com/contain-rs/bit-vec](https://github.com/contain-rs/bit-vec) |
| bitflags | 2.13.2 | MIT OR Apache-2.0 | [github.com/bitflags/bitflags](https://github.com/bitflags/bitflags) |
| block-buffer | 0.12.1 | MIT OR Apache-2.0 | [github.com/RustCrypto/utils](https://github.com/RustCrypto/utils) |
| bstr | 1.13.1 | MIT OR Apache-2.0 | [github.com/BurntSushi/bstr](https://github.com/BurntSushi/bstr) |
| bytemuck | 1.25.2 | Zlib OR Apache-2.0 OR MIT | [github.com/Lokathor/bytemuck](https://github.com/Lokathor/bytemuck) |
| bytemuck_derive | 1.12.1 | Zlib OR Apache-2.0 OR MIT | [github.com/Lokathor/bytemuck](https://github.com/Lokathor/bytemuck) |
| byteorder-lite | 0.1.0 | Unlicense OR MIT | [github.com/image-rs/byteorder-lite](https://github.com/image-rs/byteorder-lite) |
| bytes | 1.12.1 | MIT | [github.com/tokio-rs/bytes](https://github.com/tokio-rs/bytes) |
| bzip2 | 0.6.1 | MIT OR Apache-2.0 | [github.com/trifectatechfoundation/bzip2-rs](https://github.com/trifectatechfoundation/bzip2-rs) |
| cfg-if | 1.0.5 | MIT OR Apache-2.0 | [github.com/rust-lang/cfg-if](https://github.com/rust-lang/cfg-if) |
| clap | 4.6.7 | MIT OR Apache-2.0 | [github.com/clap-rs/clap](https://github.com/clap-rs/clap) |
| clap_builder | 4.6.7 | MIT OR Apache-2.0 | [github.com/clap-rs/clap](https://github.com/clap-rs/clap) |
| clap_derive | 4.6.7 | MIT OR Apache-2.0 | [github.com/clap-rs/clap](https://github.com/clap-rs/clap) |
| clap_lex | 1.1.1 | MIT OR Apache-2.0 | [github.com/clap-rs/clap](https://github.com/clap-rs/clap) |
| clipboard-win | 5.4.3 | BSL-1.0 | [github.com/DoumanAsh/clipboard-win](https://github.com/DoumanAsh/clipboard-win) |
| codespan-reporting | 0.13.1 | Apache-2.0 | [github.com/brendanzab/codespan](https://github.com/brendanzab/codespan) |
| color | 0.3.3 | Apache-2.0 OR MIT | [github.com/linebender/color](https://github.com/linebender/color) |
| colorchoice | 1.0.5 | MIT OR Apache-2.0 | [github.com/rust-cli/anstyle.git](https://github.com/rust-cli/anstyle.git) |
| const-oid | 0.10.2 | Apache-2.0 OR MIT | [github.com/RustCrypto/formats](https://github.com/RustCrypto/formats) |
| cpufeatures | 0.3.1 | MIT OR Apache-2.0 | [github.com/RustCrypto/utils](https://github.com/RustCrypto/utils) |
| crc32fast | 1.5.4 | MIT OR Apache-2.0 | [github.com/srijs/rust-crc32fast](https://github.com/srijs/rust-crc32fast) |
| crossbeam-channel | 0.5.17 | MIT OR Apache-2.0 | [github.com/crossbeam-rs/crossbeam](https://github.com/crossbeam-rs/crossbeam) |
| crossbeam-utils | 0.8.23 | MIT OR Apache-2.0 | [github.com/crossbeam-rs/crossbeam](https://github.com/crossbeam-rs/crossbeam) |
| crypto-common | 0.2.2 | MIT OR Apache-2.0 | [github.com/RustCrypto/traits](https://github.com/RustCrypto/traits) |
| cursor-icon | 1.2.0 | MIT OR Apache-2.0 OR Zlib | [github.com/rust-windowing/cursor-icon](https://github.com/rust-windowing/cursor-icon) |
| digest | 0.11.3 | MIT OR Apache-2.0 | [github.com/RustCrypto/traits](https://github.com/RustCrypto/traits) |
| document-features | 0.2.12 | MIT OR Apache-2.0 | [github.com/slint-ui/document-features](https://github.com/slint-ui/document-features) |
| dpi | 0.1.2 | Apache-2.0 AND MIT | [github.com/rust-windowing/winit](https://github.com/rust-windowing/winit) |
| ecolor | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui](https://github.com/emilk/egui) |
| eframe | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui/tree/main/crates/eframe](https://github.com/emilk/egui/tree/main/crates/eframe) |
| egui | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui](https://github.com/emilk/egui) |
| egui-wgpu | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui/tree/main/crates/egui-wgpu](https://github.com/emilk/egui/tree/main/crates/egui-wgpu) |
| egui-winit | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui/tree/main/crates/egui-winit](https://github.com/emilk/egui/tree/main/crates/egui-winit) |
| either | 1.19.0 | MIT OR Apache-2.0 | [github.com/rayon-rs/either](https://github.com/rayon-rs/either) |
| emath | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui/tree/main/crates/emath](https://github.com/emilk/egui/tree/main/crates/emath) |
| epaint | 0.37.0 | MIT OR Apache-2.0 | [github.com/emilk/egui/tree/main/crates/epaint](https://github.com/emilk/egui/tree/main/crates/epaint) |
| epaint_default_fonts | 0.37.0 | (MIT OR Apache-2.0) AND OFL-1.1 AND Ubuntu-font-1.0 | [github.com/emilk/egui/tree/main/crates/epaint_default_fonts](https://github.com/emilk/egui/tree/main/crates/epaint_default_fonts) |
| equivalent | 1.0.3 | Apache-2.0 OR MIT | [github.com/indexmap-rs/equivalent](https://github.com/indexmap-rs/equivalent) |
| error-code | 3.4.0 | BSL-1.0 | [github.com/DoumanAsh/error-code](https://github.com/DoumanAsh/error-code) |
| euclid | 0.22.14 | MIT OR Apache-2.0 | [github.com/servo/euclid](https://github.com/servo/euclid) |
| fdeflate | 0.3.7 | MIT OR Apache-2.0 | [github.com/image-rs/fdeflate](https://github.com/image-rs/fdeflate) |
| fearless_simd | 0.4.1 | Apache-2.0 OR MIT | [github.com/linebender/fearless_simd](https://github.com/linebender/fearless_simd) |
| flate2 | 1.1.10 | MIT OR Apache-2.0 | [github.com/rust-lang/flate2-rs](https://github.com/rust-lang/flate2-rs) |
| foldhash | 0.2.0 | Zlib | [github.com/orlp/foldhash](https://github.com/orlp/foldhash) |
| font-types | 0.12.7 | MIT OR Apache-2.0 | [github.com/googlefonts/fontations](https://github.com/googlefonts/fontations) |
| glow | 0.17.0 | MIT OR Apache-2.0 OR Zlib | [github.com/grovesNL/glow](https://github.com/grovesNL/glow) |
| glutin_wgl_sys | 0.6.1 | Apache-2.0 | [github.com/rust-windowing/glutin](https://github.com/rust-windowing/glutin) |
| gpu-allocator | 0.28.0 | MIT OR Apache-2.0 | [github.com/Traverse-Research/gpu-allocator](https://github.com/Traverse-Research/gpu-allocator) |
| guillotiere | 0.7.0 | MIT/Apache-2.0 | [github.com/nical/guillotiere](https://github.com/nical/guillotiere) |
| half | 2.7.1 | MIT OR Apache-2.0 | [github.com/VoidStarKat/half-rs](https://github.com/VoidStarKat/half-rs) |
| harfrust | 0.12.0 | MIT | [github.com/harfbuzz/harfrust](https://github.com/harfbuzz/harfrust) |
| hashbrown | 0.16.1 | MIT OR Apache-2.0 | [github.com/rust-lang/hashbrown](https://github.com/rust-lang/hashbrown) |
| hashbrown | 0.17.1 | MIT OR Apache-2.0 | [github.com/rust-lang/hashbrown](https://github.com/rust-lang/hashbrown) |
| heck | 0.5.0 | MIT OR Apache-2.0 | [github.com/withoutboats/heck](https://github.com/withoutboats/heck) |
| hybrid-array | 0.4.15 | MIT OR Apache-2.0 | [github.com/RustCrypto/hybrid-array](https://github.com/RustCrypto/hybrid-array) |
| image | 0.25.10 | MIT OR Apache-2.0 | [github.com/image-rs/image](https://github.com/image-rs/image) |
| indexmap | 2.14.2 | Apache-2.0 OR MIT | [github.com/indexmap-rs/indexmap](https://github.com/indexmap-rs/indexmap) |
| is_terminal_polyfill | 1.70.2 | MIT OR Apache-2.0 | [github.com/polyfill-rs/is_terminal_polyfill](https://github.com/polyfill-rs/is_terminal_polyfill) |
| itertools | 0.15.0 | MIT OR Apache-2.0 | [github.com/rust-itertools/itertools](https://github.com/rust-itertools/itertools) |
| itoa | 1.0.18 | MIT OR Apache-2.0 | [github.com/dtolnay/itoa](https://github.com/dtolnay/itoa) |
| khronos-egl | 6.0.0 | MIT/Apache-2.0 | [github.com/timothee-haudebourg/khronos-egl](https://github.com/timothee-haudebourg/khronos-egl) |
| kurbo | 0.13.1 | Apache-2.0 OR MIT | [github.com/linebender/kurbo](https://github.com/linebender/kurbo) |
| lexical-core | 1.0.6 | MIT/Apache-2.0 | [github.com/Alexhuszagh/rust-lexical](https://github.com/Alexhuszagh/rust-lexical) |
| lexical-parse-float | 1.0.6 | MIT/Apache-2.0 | [github.com/Alexhuszagh/rust-lexical](https://github.com/Alexhuszagh/rust-lexical) |
| lexical-parse-integer | 1.0.6 | MIT/Apache-2.0 | [github.com/Alexhuszagh/rust-lexical](https://github.com/Alexhuszagh/rust-lexical) |
| lexical-util | 1.0.7 | MIT/Apache-2.0 | [github.com/Alexhuszagh/rust-lexical](https://github.com/Alexhuszagh/rust-lexical) |
| lexical-write-float | 1.0.6 | MIT/Apache-2.0 | [github.com/Alexhuszagh/rust-lexical](https://github.com/Alexhuszagh/rust-lexical) |
| lexical-write-integer | 1.0.6 | MIT/Apache-2.0 | [github.com/Alexhuszagh/rust-lexical](https://github.com/Alexhuszagh/rust-lexical) |
| libbz2-rs-sys | 0.2.6 | bzip2-1.0.6 | [github.com/trifectatechfoundation/libbzip2-rs](https://github.com/trifectatechfoundation/libbzip2-rs) |
| libc | 0.2.191 | MIT OR Apache-2.0 | [github.com/rust-lang/libc](https://github.com/rust-lang/libc) |
| libloading | 0.8.9 | ISC | [github.com/nagisa/rust_libloading/](https://github.com/nagisa/rust_libloading/) |
| libm | 0.2.16 | MIT | [github.com/rust-lang/compiler-builtins](https://github.com/rust-lang/compiler-builtins) |
| libmimalloc-sys | 0.1.49 | MIT | [github.com/purpleprotocol/mimalloc_rust/tree/master/libmimalloc-sys](https://github.com/purpleprotocol/mimalloc_rust/tree/master/libmimalloc-sys) |
| linebender_resource_handle | 0.1.1 | Apache-2.0 OR MIT | [github.com/linebender/raw_resource_handle](https://github.com/linebender/raw_resource_handle) |
| litrs | 1.0.1 | MIT OR Apache-2.0 | [github.com/LukasKalbertodt/litrs](https://github.com/LukasKalbertodt/litrs) |
| lock_api | 0.4.14 | MIT OR Apache-2.0 | [github.com/Amanieu/parking_lot](https://github.com/Amanieu/parking_lot) |
| log | 0.4.35 | MIT OR Apache-2.0 | [github.com/rust-lang/log](https://github.com/rust-lang/log) |
| lzma-rust2 | 0.22.0 | Apache-2.0 | [github.com/hasenbanck/lzma-rust2/](https://github.com/hasenbanck/lzma-rust2/) |
| md-5 | 0.11.0 | MIT OR Apache-2.0 | [github.com/RustCrypto/hashes](https://github.com/RustCrypto/hashes) |
| memchr | 2.8.3 | Unlicense OR MIT | [github.com/BurntSushi/memchr](https://github.com/BurntSushi/memchr) |
| mimalloc | 0.1.52 | MIT | [github.com/purpleprotocol/mimalloc_rust](https://github.com/purpleprotocol/mimalloc_rust) |
| miniz_oxide | 0.8.9 | MIT OR Zlib OR Apache-2.0 | [github.com/Frommi/miniz_oxide/tree/master/miniz_oxide](https://github.com/Frommi/miniz_oxide/tree/master/miniz_oxide) |
| miniz_oxide | 0.9.1 | MIT OR Zlib OR Apache-2.0 | [github.com/Frommi/miniz_oxide/tree/master/miniz_oxide](https://github.com/Frommi/miniz_oxide/tree/master/miniz_oxide) |
| moxcms | 0.8.1 | BSD-3-Clause OR Apache-2.0 | [github.com/awxkee/moxcms.git](https://github.com/awxkee/moxcms.git) |
| naga | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| naga-types | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| nohash-hasher | 0.2.0 | Apache-2.0 OR MIT | [github.com/paritytech/nohash-hasher](https://github.com/paritytech/nohash-hasher) |
| num-traits | 0.2.19 | MIT OR Apache-2.0 | [github.com/rust-num/num-traits](https://github.com/rust-num/num-traits) |
| once_cell | 1.21.4 | MIT OR Apache-2.0 | [github.com/matklad/once_cell](https://github.com/matklad/once_cell) |
| once_cell_polyfill | 1.70.2 | MIT OR Apache-2.0 | [github.com/polyfill-rs/once_cell_polyfill](https://github.com/polyfill-rs/once_cell_polyfill) |
| ordered-float | 5.5.1 | MIT | [github.com/reem/rust-ordered-float](https://github.com/reem/rust-ordered-float) |
| parking_lot | 0.12.5 | MIT OR Apache-2.0 | [github.com/Amanieu/parking_lot](https://github.com/Amanieu/parking_lot) |
| parking_lot_core | 0.9.12 | MIT OR Apache-2.0 | [github.com/Amanieu/parking_lot](https://github.com/Amanieu/parking_lot) |
| pasta | 0.2.0 | MIT | [github.com/DylanCosto/pasta](https://github.com/DylanCosto/pasta) |
| peniko | 0.6.1 | Apache-2.0 OR MIT | [github.com/linebender/peniko](https://github.com/linebender/peniko) |
| pin-project-lite | 0.2.17 | Apache-2.0 OR MIT | [github.com/taiki-e/pin-project-lite](https://github.com/taiki-e/pin-project-lite) |
| png | 0.18.1 | MIT OR Apache-2.0 | [github.com/image-rs/image-png](https://github.com/image-rs/image-png) |
| pollster | 1.0.1 | Apache-2.0/MIT | [github.com/zesterer/pollster](https://github.com/zesterer/pollster) |
| polycool | 0.4.0 | MIT OR Apache-2.0 | [github.com/linebender/kurbo](https://github.com/linebender/kurbo) |
| presser | 0.3.1 | MIT OR Apache-2.0 | [github.com/EmbarkStudios/presser](https://github.com/EmbarkStudios/presser) |
| proc-macro2 | 1.0.107 | MIT OR Apache-2.0 | [github.com/dtolnay/proc-macro2](https://github.com/dtolnay/proc-macro2) |
| profiling | 1.0.18 | MIT OR Apache-2.0 | [github.com/aclysma/profiling](https://github.com/aclysma/profiling) |
| pxfm | 0.1.30 | BSD-3-Clause OR Apache-2.0 | [github.com/awxkee/pxfm](https://github.com/awxkee/pxfm) |
| quote | 1.0.47 | MIT OR Apache-2.0 | [github.com/dtolnay/quote](https://github.com/dtolnay/quote) |
| range-alloc | 0.1.5 | MIT OR Apache-2.0 | [github.com/gfx-rs/range-alloc](https://github.com/gfx-rs/range-alloc) |
| raw-window-handle | 0.6.2 | MIT OR Apache-2.0 OR Zlib | [github.com/rust-windowing/raw-window-handle](https://github.com/rust-windowing/raw-window-handle) |
| read-fonts | 0.41.0 | MIT OR Apache-2.0 | [github.com/googlefonts/fontations](https://github.com/googlefonts/fontations) |
| renderdoc-sys | 1.1.0 | MIT OR Apache-2.0 | [github.com/ebkalderon/renderdoc-rs](https://github.com/ebkalderon/renderdoc-rs) |
| rfd | 0.17.2 | MIT | [github.com/PolyMeilex/rfd](https://github.com/PolyMeilex/rfd) |
| rustc-hash | 1.1.0 | Apache-2.0/MIT | [github.com/rust-lang-nursery/rustc-hash](https://github.com/rust-lang-nursery/rustc-hash) |
| rustc-hash | 2.1.3 | Apache-2.0 OR MIT | [github.com/rust-lang/rustc-hash](https://github.com/rust-lang/rustc-hash) |
| scopeguard | 1.2.0 | MIT OR Apache-2.0 | [github.com/bluss/scopeguard](https://github.com/bluss/scopeguard) |
| self_cell | 1.3.0 | Apache-2.0 OR GPL-2.0-only | [github.com/Voultapher/self_cell](https://github.com/Voultapher/self_cell) |
| serde | 1.0.229 | MIT OR Apache-2.0 | [github.com/serde-rs/serde](https://github.com/serde-rs/serde) |
| serde_core | 1.0.229 | MIT OR Apache-2.0 | [github.com/serde-rs/serde](https://github.com/serde-rs/serde) |
| serde_derive | 1.0.229 | MIT OR Apache-2.0 | [github.com/serde-rs/serde](https://github.com/serde-rs/serde) |
| serde_json | 1.0.154 | MIT OR Apache-2.0 | [github.com/serde-rs/json](https://github.com/serde-rs/json) |
| sha2 | 0.11.1 | MIT OR Apache-2.0 | [github.com/RustCrypto/hashes](https://github.com/RustCrypto/hashes) |
| simd-adler32 | 0.3.10 | MIT | [github.com/mcountryman/simd-adler32](https://github.com/mcountryman/simd-adler32) |
| skrifa | 0.44.0 | MIT OR Apache-2.0 | [github.com/googlefonts/fontations](https://github.com/googlefonts/fontations) |
| smallvec | 1.16.4 | MIT OR Apache-2.0 | [github.com/servo/rust-smallvec](https://github.com/servo/rust-smallvec) |
| smol_str | 0.2.2 | MIT OR Apache-2.0 | [github.com/rust-analyzer/smol_str](https://github.com/rust-analyzer/smol_str) |
| spirv | 0.4.0+sdk-1.4.341.0 | Apache-2.0 | [github.com/gfx-rs/rspirv](https://github.com/gfx-rs/rspirv) |
| static_assertions | 1.1.0 | MIT OR Apache-2.0 | [github.com/nvzqz/static-assertions-rs](https://github.com/nvzqz/static-assertions-rs) |
| strict-num | 0.1.1 | MIT | [github.com/RazrFalcon/strict-num](https://github.com/RazrFalcon/strict-num) |
| strsim | 0.11.1 | MIT | [github.com/rapidfuzz/strsim-rs](https://github.com/rapidfuzz/strsim-rs) |
| syn | 2.0.119 | MIT OR Apache-2.0 | [github.com/dtolnay/syn](https://github.com/dtolnay/syn) |
| syn | 3.0.9 | MIT OR Apache-2.0 | [github.com/dtolnay/syn](https://github.com/dtolnay/syn) |
| thiserror | 2.0.21 | MIT OR Apache-2.0 | [github.com/dtolnay/thiserror](https://github.com/dtolnay/thiserror) |
| thiserror-impl | 2.0.21 | MIT OR Apache-2.0 | [github.com/dtolnay/thiserror](https://github.com/dtolnay/thiserror) |
| tiny-skia | 0.12.0 | BSD-3-Clause | [github.com/linebender/tiny-skia](https://github.com/linebender/tiny-skia) |
| tiny-skia-path | 0.12.0 | BSD-3-Clause | [github.com/linebender/tiny-skia/tree/master/path](https://github.com/linebender/tiny-skia/tree/master/path) |
| tracing | 0.1.44 | MIT | [github.com/tokio-rs/tracing](https://github.com/tokio-rs/tracing) |
| tracing-core | 0.1.36 | MIT | [github.com/tokio-rs/tracing](https://github.com/tokio-rs/tracing) |
| type-map | 0.5.1 | MIT/Apache-2.0 | [github.com/kardeiz/type-map](https://github.com/kardeiz/type-map) |
| typenum | 1.20.1 | MIT OR Apache-2.0 | [github.com/paholg/typenum](https://github.com/paholg/typenum) |
| unicode-general-category | 1.1.0 | Apache-2.0 | [github.com/yeslogic/unicode-general-category](https://github.com/yeslogic/unicode-general-category) |
| unicode-ident | 1.0.27 | (MIT OR Apache-2.0) AND Unicode-3.0 | [github.com/dtolnay/unicode-ident](https://github.com/dtolnay/unicode-ident) |
| unicode-properties | 0.1.5 | (MIT OR Apache-2.0) AND Unicode-3.0 | [github.com/unicode-rs/unicode-properties](https://github.com/unicode-rs/unicode-properties) |
| unicode-segmentation | 1.13.3 | MIT OR Apache-2.0 | [github.com/unicode-rs/unicode-segmentation](https://github.com/unicode-rs/unicode-segmentation) |
| unicode-width | 0.2.3 | MIT OR Apache-2.0 | [github.com/unicode-rs/unicode-width](https://github.com/unicode-rs/unicode-width) |
| utf8parse | 0.2.2 | Apache-2.0 OR MIT | [github.com/alacritty/vte](https://github.com/alacritty/vte) |
| uuid | 1.29.0 | Apache-2.0 OR MIT | [github.com/uuid-rs/uuid](https://github.com/uuid-rs/uuid) |
| vello_common | 0.2.0 | Apache-2.0 OR MIT | [github.com/linebender/vello](https://github.com/linebender/vello) |
| vello_cpu | 0.2.0 | Apache-2.0 OR MIT | [github.com/linebender/vello](https://github.com/linebender/vello) |
| web-time | 1.1.0 | MIT OR Apache-2.0 | [github.com/daxpedda/web-time](https://github.com/daxpedda/web-time) |
| wgpu | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| wgpu-core | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| wgpu-core-deps-windows-linux-android | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| wgpu-hal | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| wgpu-naga-bridge | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| wgpu-types | 30.0.1 | MIT OR Apache-2.0 | [github.com/gfx-rs/wgpu](https://github.com/gfx-rs/wgpu) |
| windows | 0.62.2 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-collections | 0.3.2 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-core | 0.62.2 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-future | 0.3.2 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-implement | 0.60.2 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-interface | 0.59.3 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-link | 0.2.1 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-numerics | 0.3.1 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-result | 0.4.1 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-strings | 0.5.1 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-sys | 0.52.0 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-sys | 0.61.2 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-targets | 0.52.6 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows-threading | 0.2.1 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| windows_x86_64_msvc | 0.52.6 | MIT OR Apache-2.0 | [github.com/microsoft/windows-rs](https://github.com/microsoft/windows-rs) |
| winit | 0.30.13 | Apache-2.0 | [github.com/rust-windowing/winit](https://github.com/rust-windowing/winit) |
| zerocopy | 0.8.63 | BSD-2-Clause OR Apache-2.0 OR MIT | [github.com/google/zerocopy](https://github.com/google/zerocopy) |
| zerocopy-derive | 0.8.63 | BSD-2-Clause OR Apache-2.0 OR MIT | [github.com/google/zerocopy](https://github.com/google/zerocopy) |
| zlib-rs | 0.6.8 | Zlib | [github.com/trifectatechfoundation/zlib-rs](https://github.com/trifectatechfoundation/zlib-rs) |
| zmij | 1.0.23 | MIT | [github.com/dtolnay/zmij](https://github.com/dtolnay/zmij) |
