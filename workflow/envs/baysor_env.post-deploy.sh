
# workaround for missing $CONDA_PREFIX/share/julia/cert.pem
# see https://bugs.gentoo.org/888978
if [ ! -e "$CONDA_PREFIX/share/julia/cert.pem" ] ; then
	ln -s "$CONDA_PREFIX/ssl/cert.pem" "$CONDA_PREFIX/share/julia/cert.pem"
fi

# see https://discourse.julialang.org/t/using-plots-fails-partially-solved-but-mysterious/76120
echo "" >> "$CONDA_PREFIX/etc/conda/activate.d/julia_activate.sh"
echo "# workaround for incompatible LD_LIBRARY_PATH settings" >> "$CONDA_PREFIX/etc/conda/activate.d/julia_activate.sh"
echo "alias julia='LD_LIBRARY_PATH=\"\" \"$CONDA_PREFIX/bin/julia\" -q'" >> "$CONDA_PREFIX/etc/conda/activate.d/julia_activate.sh"

# install Baysor
LD_LIBRARY_PATH="" "$CONDA_PREFIX/bin/julia" -q -e 'using Pkg; Pkg.add([PackageSpec(name="UMAP", rev="master"), PackageSpec(url="https://github.com/kharchenkolab/Baysor.git")]); Pkg.add(["Colors","Plots","CSV","Distances","ProgressMeter","DataFrames","DataFramesMeta","NearestNeighbors","RCall","StatsBase","ArgParse"]); import Baysor'
