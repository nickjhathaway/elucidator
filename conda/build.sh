#!/usr/bin/env bash
# Builds elucidator and its dependencies with setup.py into $PREFIX/share/elucidator,
# trims what isn't needed at runtime, and links the elucidator, SeekDeep and PathWeaver
# executables into $PREFIX/bin.
set -euxo pipefail

if [[ "${target_platform}" == osx-* ]]; then
  SETUP_CC="${CC}"
  SETUP_CXX="${CXX}"
  # restbed's bundled asio needs <exception>
  export CXXFLAGS="${CXXFLAGS:-} -include exception"
else
  # setup.py expects compilers named gcc-NN / g++-NN
  GCCVER=$("${CXX}" -dumpversion | cut -d. -f1)
  WRAP="${SRC_DIR}/_toolwrap"
  mkdir -p "${WRAP}"
  printf '#!/usr/bin/env bash\nexec %s "$@"\n' "$(command -v "${CC}")" > "${WRAP}/gcc-${GCCVER}"
  printf '#!/usr/bin/env bash\nexec %s "$@"\n' "$(command -v "${CXX}")" > "${WRAP}/g++-${GCCVER}"
  chmod +x "${WRAP}/gcc-${GCCVER}" "${WRAP}/g++-${GCCVER}"
  export PATH="${WRAP}:${PATH}"
  SETUP_CC="gcc-${GCCVER}"
  SETUP_CXX="g++-${GCCVER}"
fi

INST="${PREFIX}/share/elucidator"
mkdir -p "${INST}"
cd "${INST}"
python "${SRC_DIR}/setup.py" \
  --CC "${SETUP_CC}" --CXX "${SETUP_CXX}" \
  --libs "elucidator:develop:${GIT_REV:-82a4f0f}" \
  --numCores "${CPU_COUNT:-8}" --overWrite

# keep only what's needed at runtime
rm -rf external/build external/tarballs external/.cache
find external/local \( -name .git -o -name tests -o -name docs -o -name doc \) -type d -prune -exec rm -rf {} +
find external/local -name '*.o' -delete
find external/local -mindepth 4 -maxdepth 4 -type d -name include -prune -exec rm -rf {} +
find external/local \( -name '*.a' -o -name '*.la' \) -type f -delete
find external/local \( -path '*/lib/cmake' -o -path '*/lib/pkgconfig' -o -path '*/share/man' \) -type d -prune -exec rm -rf {} +
find external/local -depth -type d -empty -delete

# fail if any binary needs a newer glibc / macOS than TARGET_STDLIB_VERSION
if [ -n "${TARGET_STDLIB_VERSION:-}" ]; then
  if [[ "${target_platform}" == osx-* ]]; then
    what="macOS"
    needed=$(find external/local -type f \( -name '*.dylib' -o -perm -u+x \) \
               -exec "${OTOOL:-otool}" -l {} + 2>/dev/null \
             | awk '/^ *minos / {print $2}' | sort -uV | tail -1 || true)
  else
    what="glibc"
    needed=$(find external/local -type f \( -name '*.so*' -o -perm -u+x \) \
               -exec "${OBJDUMP:-objdump}" -T {} + 2>/dev/null \
             | grep -oE 'GLIBC_[0-9]+(\.[0-9]+)+' | sed 's/^GLIBC_//' | sort -uV | tail -1 || true)
  fi
  echo "highest ${what} version needed: ${needed:-none} (target ${TARGET_STDLIB_VERSION})"
  if [ -n "${needed}" ] && \
     [ "$(printf '%s\n%s\n' "${needed}" "${TARGET_STDLIB_VERSION}" | sort -V | tail -1)" != "${TARGET_STDLIB_VERSION}" ]; then
    echo "ERROR: binaries need ${what} ${needed}, newer than the ${TARGET_STDLIB_VERSION} target" >&2
    exit 1
  fi
fi

# put only the njh tools on PATH
mkdir -p "${PREFIX}/bin"
for tool in elucidator SeekDeep PathWeaver; do
  for exe in external/local/"${tool}"/*/"${tool}"/bin/*; do
    [ -f "${exe}" ] && [ -x "${exe}" ] || continue
    ln -sf "../share/elucidator/${exe}" "${PREFIX}/bin/$(basename "${exe}")"
  done
done
