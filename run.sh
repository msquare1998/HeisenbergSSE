lx=16
ly=1      # set this to be "1" for simulating 1D chain
beta=16
nThm=100000
nStat=50000
nBins=10
if [[ "$(uname)" == "Darwin" ]]; then
    echo "- No extra option selected on macOS."
elif [[ "$(uname)" == "Linux" ]]; then
    RUSTFLAGS="-C target-feature=+avx,+fma" cargo build --release
fi
cargo build --release
./target/release/HM_2D $lx $ly $beta $nThm $nStat $nBins