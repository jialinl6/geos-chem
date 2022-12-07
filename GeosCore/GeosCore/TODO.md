# TODO

## Initial transliteration

This involves, initial translation, and documentation. 100% match in functionality minus CPU flags & macros & compiler flags.

- [ ] `./src/tpcore_fvdas_mod.jl`
  - [ ] `init_tpcore!`
  - [ ] `tpcore_fvdas`
  - [x] `average_const_poles!`
  - [x] `set_cross_terms!`
  - [x] `calc_vert_mass_flux!`
  - [x] `set_jn_js!`
  - [x] `calc_advec_cross_terms!`
  - [x] `qckxyz!`
  - [x] `set_lmts!`
  - [x] `set_press_terms!`
  - [x] `calc_courant!`
  - [x] `calc_divergence!`
  - [x] `do_divergence_pole_sum!`
  - [x] `do_cross_terms_pole_i2d2!`
  - [x] `xadv_dao2!`
  - [x] `yadv_dao2!`
  - [x] `do_yadv_pole_i2d2!`
  - [x] `do_yadv_pole_sum!`
  - [x] `xtp!`
  - [x] `xmist!`
  - [x] `fxppm!`
  - [x] `lmtppm!`
  - [x] `ytp!`
  - [x] `ymist!`
  - [x] `do_ymist_pole1_i2d2!`
  - [x] `do_ymist_pole2_i2d2!`
  - [x] `fyppm!`
  - [x] `do_fyppm_pole_i2d2!`
  - [x] `do_ytp_pole_sum!`
  - [x] `fzppm!`
  - [x] `average_press_poles!`
- [ ] `./src/tpcore_window_mod.jl`
- [ ] `./src/transport_mod.jl`
- [ ] `./src/drydep_mod.jl`
- [ ] `./src/aerosol_mod.jl`
