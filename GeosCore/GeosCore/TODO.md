# TODO

## Initial transliteration

This involves, initial translation, and documentation. 100% match in functionality minus CPU flags & macros & compiler flags.

- [ ] `./tpcore_fvdas_mod.jl`
  - [ ] `init_tpcore!`
  - [ ] `tpcore_fvdas`
  - [ ] `average_const_poles!`
  - [ ] `set_cross_terms!`
  - [ ] `calc_vert_mass_flux!`
  - [ ] `set_jn_js!`
  - [ ] `calc_advec_cross_terms!`
  - [ ] `qckxyz!`
  - [ ] `set_lmts!`
  - [ ] `set_press_terms!`
  - [ ] `calc_courant!`
  - [ ] `calc_divergence!`
  - [ ] `do_divergence_pole_sum!`
  - [ ] `do_cross_terms_pole_i2d2!`
  - [ ] `xadv_dao2!`
  - [ ] `yadv_dao2!`
  - [ ] `do_yadv_pole_i2d2!`
  - [ ] `do_yadv_pole_sum!`
  - [ ] `xtp!`
  - [ ] `xmist!`
  - [ ] `fxppm!`
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
- [ ] `./tpcore_window_mod.jl`
- [ ] `./transport_mod.jl`
- [ ] `./drydep_mod.jl`
- [ ] `./aerosol_mod.jl`
