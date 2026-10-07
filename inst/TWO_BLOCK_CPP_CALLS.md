# Two-block C++ symbols and the R migration

The groups/population rename is applied in glmbayesCore C++. lmerb calls three
of these exports from R: `two_block_build_mu_all_cpp_export`,
`two_block_population_one_chain_cpp_export`, and
`rIndepNormalGammaRegBlockInd_cpp_export`. The others are marked
`lmerb-unused` in the `.cpp` sources and are candidates for deletion.
lmebayesCore wrappers still use the old names; do not copy them again.

Argument order is unchanged. Argument names `prior_list_block1`,
`dispersion_block1`, `ddef_block1`, and `block1_prior` are unchanged.

## Retarget these `.Call` symbols

| lmebayesCore R wrapper | lmebayesCore `.Call` symbol | glmbayesCore `.Call` symbol |
|---|---|---|
| `.two_block_rNormal_reg_cpp` | `two_block_rNormal_reg_v2_cpp_export` | `two_block_rNormal_reg_cpp_export` |
| `.two_block_block1_prior_with_tau2_cpp` | `two_block_block1_prior_with_tau2_cpp_export` | `two_block_groups_prior_with_tau2_cpp_export` |
| `.two_block_block1_iters_mean_cpp` | `two_block_block1_iters_mean_cpp_export` | `two_block_groups_iters_mean_cpp_export` |
| `.two_block_block1_one_chain_draw_cpp` | `two_block_block1_one_chain_draw_cpp_export` | `two_block_groups_one_chain_draw_cpp_export` |
| `.two_block_block1_one_chain_cpp` | `two_block_block1_one_chain_cpp_export` | `two_block_groups_one_chain_cpp_export` |
| `.two_block_block1_one_chain_from_mu_P_cpp` | `two_block_block1_one_chain_from_mu_P_cpp_export` | `two_block_groups_one_chain_from_mu_P_cpp_export` |
| `.two_block_block2_one_chain_cpp` | `two_block_block2_one_chain_cpp_export` | `two_block_population_one_chain_cpp_export` |

`.two_block_rNormal_reg_cpp` must call
`_glmbayesCore_two_block_rNormal_reg_cpp_export`. The `_v2` suffix was
dropped because glmbayesCore has no other `two_block_rNormal_reg`. Do not
copy lmebayesCore's generated `two_block_rNormal_reg_cpp_export`: that symbol
there is the removed v1 driver and takes `prior_list_block2`. The
glmbayesCore symbol of the same name takes `pfamily_list`.

## Keep, and change only the package prefix

| lmebayesCore R wrapper | `.Call` symbol in both packages |
|---|---|
| `.two_block_build_mu_all_cpp` | `two_block_build_mu_all_cpp_export` |
| `.two_block_batch_fixef_chain_cpp` | `two_block_batch_fixef_chain_cpp_export` |
| `.two_block_batch_tau2_chain_row_cpp` | `two_block_batch_tau2_chain_row_cpp_export` |
| `.two_block_batch_b_assign_slice_cpp` | `two_block_batch_b_assign_slice_cpp_export` |
| `.two_block_batch_iters_ranef_add_cpp` | `two_block_batch_iters_ranef_add_cpp_export` |
| `.two_block_reorder_b_to_group_levels_cpp` | `two_block_reorder_b_to_group_levels_cpp_export` |
| `.two_block_align_b_to_xhyper_cpp` | `two_block_align_b_to_xhyper_cpp_export` |
| `.BlockEnvelopeCentering_cpp` | `BlockEnvelopeCentering_cpp_export` |
| `.BlockEnvelopeBuild_cpp` | `BlockEnvelopeBuild_cpp_export` |
| `.BlockEnvelopeDispersionBuild_cpp` | `BlockEnvelopeDispersionBuild_cpp_export` |
| `.BlockEnvelopeSim_cpp` | `BlockEnvelopeSim_cpp_export` |
| `.rIndepNormalGammaRegBlock_cpp` | `rIndepNormalGammaRegBlock_cpp_export` |
| `.BlockEnvelopeDispersionBuildInd_cpp` | `BlockEnvelopeDispersionBuildInd_cpp_export` |
| `.BlockEnvelopeSimInd_cpp` | `BlockEnvelopeSimInd_cpp_export` |
| `.rIndepNormalGammaRegBlockInd_cpp` | `rIndepNormalGammaRegBlockInd_cpp_export` |

## Do not copy these call sites

These symbols were removed from glmbayesCore. They remain in lmebayesCore.
The R functions that called them have no caller on the current sweep path.

| lmebayesCore R wrapper | `.Call` symbol | Note |
|---|---|---|
| generated v1 `two_block_rNormal_reg_cpp_export` | same spelling, different function | lmebayesCore's unsuffixed export is the v1 driver (`prior_list_block2`). glmbayesCore reused that spelling for the pfamily driver. Copy `.two_block_rNormal_reg_cpp`, not the generated v1 wrapper. |
| `.two_block_rNormal_reg_staged_cpp` | `two_block_rNormal_reg_staged_cpp_export` | Staged driver was not ported. |
| generated `two_block_rNormal_reg_v5_cpp_export` | `two_block_rNormal_reg_v5_cpp_export` | v5 sweep-outer driver. v3 and v4 were never exported. |
| `.two_block_block1_one_chain_v2_cpp` | `two_block_block1_one_chain_v2_cpp_export` | Unused groups-update shape. |
| `.two_block_block1_all_chains_v2_internal_cpp` | `two_block_block1_all_chains_v2_internal_cpp_export` | Unused batch loop. |
| `.two_block_block1_all_chains_v2_internal_loop_cpp` | `two_block_block1_all_chains_v2_internal_loop_cpp_export` | Unused batch loop. |
| `.two_block_block1_all_chains_cpp` | `two_block_block1_all_chains_cpp_export` | Unused batch loop. Called only from `.two_block_block1_all_chains_via_cpp`, which itself has no caller. |

## Internal C++ names, no R `.Call`

| Former C++ name | Current C++ name |
|---|---|
| `block1_prior_list` | `groups_prior_list` |
| `two_block_block2_rglmb_gamma` | `two_block_population_rglmb_gamma` |
| `block2_prior_prep_v2` | `population_prior_prep` |
| `block2_mu_P_prep_v2` | `population_mu_P_prep` |
| `Block2PriorV2` | `PopulationPrior` |

`src/two_block_block1.cpp` is now `src/two_block_groups.cpp`. The sweep
population helpers are in `src/two_block_population.cpp`. The Gibbs driver
stays in `src/twoBlockGibbs.cpp`. Those three files are in
`glmbayes::sim::twoblock`. `block_rIndepNormalGammaReg.cpp` is in
`glmbayes::env::group`. The `.Call` symbol names are unchanged by that move.

These names do not contain `block1` or `block2` and stay:
`two_block_rNormal_reg`, `two_block_build_mu_all`, `two_block_batch_*`,
`two_block_reorder_b_to_group_levels`, `two_block_align_b_to_xhyper`, and the
`BlockEnvelope*` / `rIndepNormalGammaRegBlock*` entry points.
