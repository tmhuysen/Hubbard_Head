# This file specifies PACKAGE_VERSION for the "Hubbard_Head" package.

# It sets
#   - PACKAGE_VERSION_EXACT         to TRUE     if the requested version is equal to the current version
#   - PACKAGE_VERSION_COMPATIBLE    to FALSE    if the current version < therequested version
#   - PACKAGE_VERSION_COMPATIBLE    to TRUE     if the current version >= the requested version


set(PACKAGE_VERSION 0.0.1)

if("${PACKAGE_VERSION}" VERSION_LESS "${PACKAGE_FIND_VERSION}")
    set(PACKAGE_VERSION_COMPATIBLE FALSE)
else()
    set(PACKAGE_VERSION_COMPATIBLE TRUE)
    if("${PACKAGE_FIND_VERSION}" STREQUAL "${PACKAGE_VERSION}")
        set(PACKAGE_VERSION_EXACT TRUE)
    endif()
endif()
