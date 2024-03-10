#curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
export MAMBA_ROOT_PREFIX=$SCRATCH/micromamba # optional, defaults to ~/micromamba
pushd $MAMBA_ROOT_PREFIX
  eval "$(./bin/micromamba shell hook -s posix)"
  micromamba activate
popd
#micromamba install -c conda-forge root=6.18
