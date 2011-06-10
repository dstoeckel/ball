### list all filenames of the directory here ###
SET(SOURCES_LIST
argedit.C
argloader.C
argraph.C
error.C
gene.C
gene_mesh.C
match.C
sd_state.C
sortnodes.C
ull_state.C
ull_sub_state.C
vf2_mcs_state.C
vf2_mono_state.C
vf2_state.C
vf2_sub_state.C
vf_mono_state.C
vf_state.C
vf_sub_state.C
xsubgraph.C
)

ADD_BALL_SOURCES("STRUCTURE/VFLIB" "${SOURCES_LIST}")
