#!/usr/bin/env bash
set -euo pipefail

# -------------------------------------------------------
# Handle --clean option FIRST
# -------------------------------------------------------
if [[ "${1:-}" == "--clean" ]]; then
  ROOT_DIR="$(pwd)"

  if [[ "${2:-}" == "site" ]]; then
    echo "Cleaning build/site ..."
    rm -rf "${ROOT_DIR}/build/site"
  elif [[ "${2:-}" == "global" ]]; then
    echo "Cleaning build/global ..."
    rm -rf "${ROOT_DIR}/build/global"
  else
    echo "Cleaning entire build directory ..."
    rm -rf "${ROOT_DIR}/build"
  fi

  echo "Clean complete."
  exit 0
fi

# -------------------------------------------------------
# Defaults
# -------------------------------------------------------
DEFAULT_MODE="site"
DEFAULT_NML="./para_files/parameters_ORNL_test.nml"

#DEFAULT_MODE="global"
#DEFAULT_NML="./para_files/parameters_Global_test.nml"

CPPFLAGS=( -DScreenOutput )
LDFLAGS=()

# -------------------------------------------------------
# NetCDF default locations (edit for your machine)
# Can still be overridden by environment variables
# -------------------------------------------------------
#NETCDF_INC="${NETCDF_INC:-}"
#NETCDF_LIB="${NETCDF_LIB:-}"
NETCDF_INC="${NETCDF_INC:-/opt/local/include}"
NETCDF_LIB="${NETCDF_LIB:-/opt/local/lib}"

# -------------------------------------------------------
# Model mode and namelist file
MODE="${1:-$DEFAULT_MODE}"
NML="${2:-$DEFAULT_NML}"

if [[ "${MODE}" != "site" && "${MODE}" != "global" ]]; then
  echo "ERROR: MODE must be 'site' or 'global' (got: ${MODE})"
  exit 2
fi
if [[ ! -f "${NML}" ]]; then
  echo "ERROR: Namelist not found: ${NML}"
  exit 2
fi

FC="${FC:-gfortran}"

ROOT_DIR="$(pwd)"
SRC_DIR="${ROOT_DIR}/src"
BUILD_DIR="${ROOT_DIR}/build/${MODE}"
EXE="${BUILD_DIR}/ess_${MODE}"

mkdir -p "${BUILD_DIR}"
pushd "${BUILD_DIR}" >/dev/null

MODDIR="./mod"
mkdir -p "${MODDIR}"

COMMON_FFLAGS=( -O2 -cpp -ffree-line-length-none -J"${MODDIR}" -I"${MODDIR}" )
if [[ "${DEBUG:-0}" == "1" ]]; then
  COMMON_FFLAGS=( -O0 -g -cpp -ffree-line-length-none -fcheck=all -fbacktrace -J"${MODDIR}" -I"${MODDIR}" )
fi

# -------------------------------------------------------
# Source files (correct module order)
# -------------------------------------------------------
SRCS=(
  "${SRC_DIR}/datatypes.F90"
  "${SRC_DIR}/model_utils.F90"
  "${SRC_DIR}/io_mod.F90"
  "${SRC_DIR}/soil.F90"
  "${SRC_DIR}/vegetation.F90"
  "${SRC_DIR}/BiomeE.F90"
)

if [[ "${MODE}" == "global" ]]; then
  CPPFLAGS+=( -DGlobalRun -DDO_Climate_VEG -DZip_outputs -DZippedNCfiles -DUSE_NETCDF )
  SRCS+=( "${SRC_DIR}/netcdf_io.F90" )

  if [[ ! -f "${NETCDF_INC}/netcdf.mod" ]]; then
    echo "ERROR: netcdf.mod not found in NETCDF_INC=${NETCDF_INC}"
    echo "       Set NETCDF_INC/NETCDF_LIB or install netcdf-fortran."
    exit 2
  fi

  [[ -n "${NETCDF_INC}" ]] && COMMON_FFLAGS+=( "-I${NETCDF_INC}" )
  [[ -n "${NETCDF_LIB}" ]] && LDFLAGS+=( "-L${NETCDF_LIB}" )
  LDFLAGS+=( -lnetcdff -lnetcdf )
fi

SRCS+=( "${SRC_DIR}/main.F90" )

# -------------------------------------------------------
# Compile
# -------------------------------------------------------
echo "==> Building (${MODE})"
rm -f ./*.o "${MODDIR}"/*.mod "${MODDIR}"/*.smod "${EXE}"

CMD=( "${FC}" )
CMD+=( "${COMMON_FFLAGS[@]}" )
CMD+=( "${CPPFLAGS[@]}" )
CMD+=( "${SRCS[@]}" )
if [[ ${#LDFLAGS[@]} -gt 0 ]]; then CMD+=( "${LDFLAGS[@]}" ); fi
CMD+=( -o "${EXE}" )

"${CMD[@]}"

popd >/dev/null

# -------------------------------------------------------
# Run from project root (preserves relative paths)
# -------------------------------------------------------
mkdir -p "${ROOT_DIR}/para_files"
ln -sf "$(realpath "${NML}")" "${ROOT_DIR}/para_files/input.nml"

echo "==> Running ${EXE}"
pushd "${ROOT_DIR}" >/dev/null
"${EXE}"
popd >/dev/null
