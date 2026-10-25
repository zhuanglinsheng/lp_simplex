# PowerShell 7.3+; dependencies live in project-owned build directories.
param([string]$VcpkgRoot = (Join-Path $PSScriptRoot '../build-vcpkg'))
$ErrorActionPreference = 'Stop'
$PSNativeCommandUseErrorActionPreference = $true
$root = (Resolve-Path (Join-Path $PSScriptRoot '..')).Path
$VcpkgRoot = [IO.Path]::GetFullPath($VcpkgRoot)
$manifest = Get-Content (Join-Path $root 'vcpkg.json') -Raw | ConvertFrom-Json
$revision = $manifest.'builtin-baseline'
if (!(Test-Path (Join-Path $VcpkgRoot '.git'))) {
    if (Test-Path $VcpkgRoot) {
        throw "Refusing to replace existing directory: $VcpkgRoot"
    }
    git init $VcpkgRoot
    git -C $VcpkgRoot remote add origin https://github.com/microsoft/vcpkg.git
    git -C $VcpkgRoot fetch --depth 1 origin $revision
    git -C $VcpkgRoot checkout --detach $revision
} elseif ((git -C $VcpkgRoot rev-parse HEAD).Trim() -ne $revision) {
    throw 'Existing vcpkg checkout has another revision; choose a new VcpkgRoot.'
}
& (Join-Path $VcpkgRoot 'bootstrap-vcpkg.bat') -disableMetrics
& (Join-Path $VcpkgRoot 'vcpkg.exe') install --triplet x64-windows `
    "--x-manifest-root=$root" "--x-install-root=$root/build-vcpkg-installed"
