# Run from a Visual Studio x64 developer shell with PowerShell 7.3+.
param(
    [string]$BuildDir = (Join-Path $PSScriptRoot '../build-windows'),
    [string]$VcpkgRoot = (Join-Path $PSScriptRoot '../build-vcpkg'),
    [ValidateSet('Release', 'Debug')][string]$Configuration = 'Release',
    [switch]$Minimal
)
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$BuildDir = [IO.Path]::GetFullPath($BuildDir)
$VcpkgRoot = [IO.Path]::GetFullPath($VcpkgRoot)
$toolchain = Join-Path $VcpkgRoot 'scripts/buildsystems/vcpkg.cmake'
$installed = Join-Path $root 'build-vcpkg-installed'
if (!(Test-Path $toolchain)) {
    throw 'Run tools/install_dependencies.ps1 first.'
}
$backends = if ($Minimal) { 'OFF' } else { 'ON' }
cmake -S $root -B $BuildDir -G Ninja "-DCMAKE_BUILD_TYPE=$Configuration" `
    "-DCMAKE_TOOLCHAIN_FILE=$toolchain" "-DVCPKG_INSTALLED_DIR=$installed" `
    -DVCPKG_TARGET_TRIPLET=x64-windows -DLP_SIMPLEX_REQUIRE_OPTIONAL_BACKENDS=ON `
    "-DLP_SIMPLEX_USE_KLU=$backends" "-DLP_SIMPLEX_USE_SPQR=$backends" `
    "-DLP_SIMPLEX_USE_GMP=$backends"
cmake --build $BuildDir --parallel 3
$env:PATH = "$BuildDir;$installed/x64-windows/bin;$env:PATH"
if ($Configuration -eq 'Debug') {
    $env:PATH = "$installed/x64-windows/debug/bin;$env:PATH"
}
ctest --test-dir $BuildDir --output-on-failure
