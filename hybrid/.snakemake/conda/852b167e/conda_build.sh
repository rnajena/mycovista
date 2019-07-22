if [ -z ${CONDA_BUILD+x} ]; then
	source /opt/conda/conda-bld/circos_1562082357450/work/build_env_setup.sh
fi
#!/bin/bash

cp -r . $PREFIX
