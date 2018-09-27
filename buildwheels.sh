#!/bin/bash
#
# Build manylinux1 wheels. Based on the example at
# <https://github.com/pypa/python-manylinux-demo>
#
# For interactive tests:
#   docker run -it -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /bin/bash

set -xeuo pipefail

# For convenience, if this script is called from outside of a docker container,
# it starts a container and runs itself inside of it.
if ! grep -q docker /proc/1/cgroup; then
  # We are not inside a container
  docker pull quay.io/pypa/manylinux1_x86_64
  exec docker run --rm -v $(pwd):/io quay.io/pypa/manylinux1_x86_64 /io/$0
fi

# We don’t support Python 2.7
rm /opt/python/cp27*

PYBINS="/opt/python/*/bin"
HAS_CYTHON=0
for PYBIN in ${PYBINS}; do
    ${PYBIN}/pip install Cython
#    ${PYBIN}/pip install -r /io/requirements.txt
    ${PYBIN}/pip wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/dnaio-*.whl; do
    auditwheel repair "$whl" -w repaired/
done

# Created files are owned by root, so fix permissions.
chown -R --reference=/io/setup.py repaired/
mv repaired/*.whl /io/dist/
