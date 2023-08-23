#!/bin/bash
# Create a Conda environment named "rust" with Python 3.7
conda create --name rust python=3.7

# Activate the "rust" Conda environment
source activate rust

# Use the Ubuntu Advanced Packaging Tool (apt) to install the build-essential package
sudo apt install build-essential

# Install various dependencies required for the following steps
sudo apt install --assume-yes git clang curl libssl-dev protobuf-compiler

# Install the Rust standard library for x86_64 architecture from the conda-forge channel
conda install -c conda-forge rust-std-x86_64-unknown-linux-gnu

# Download the rustup installation program and use it to install Rust
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

# Update your current shell to include Cargo, Rust's package manager and build tool
source $HOME/.cargo/env

# Verify your Rust installation by checking the Rust compiler version
rustc --version

# Configure the Rust toolchain to default to the latest stable version
rustup default stable

# Update the Rust toolchain to the latest stable version
rustup update

# Add the nightly release and the nightly WebAssembly (wasm) targets to your development environment
rustup update nightly
rustup target add wasm32-unknown-unknown --toolchain nightly

# Verify the configuration of your development environment
rustup show
rustup +nightly show

# Clone the "transanno" repository into the specified directory
mkdir -p DNA_softwares/transanno
cd $HOME/DNA_softwares/transanno
git clone https://github.com/informationsea/transanno.git

# Run a script called "prepare-test-files.sh" in the "transanno" directory
cd transanno
./prepare-test-files.sh

# Run tests for the "transanno" project using Cargo and build it in release mode
cargo test && cargo build --release

# Install the 'bcftools' and 'minimap2' packages from the bioconda channel, supporting liftOver
conda install -c bioconda bcftools

conda install -c bioconda minimap2

# Deactivate the "rust" Conda environment
conda deactivate

