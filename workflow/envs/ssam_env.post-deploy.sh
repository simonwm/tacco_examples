
# workaround for pyarrow
if [ ! -e "$CONDA_PREFIX/lib/libglog.so.0" ] ; then
	ln -s $CONDA_PREFIX/lib/libglog.so $CONDA_PREFIX/lib/libglog.so.0
fi
