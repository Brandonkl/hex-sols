```sh
git clone https://github.com/Brandonkl/moire-exact
git clone https://github.com/Brandonkl/hex-sols
git clone https://github.com/Brandonkl/rust-sols

(cd moire-exact && python3 -m pip install .)  # this is a library used by positions.py

# Generate all [a,b,c] triplets in (theta, 60deg - theta) pairs
# that have volume < 250
cd rust-sols
cargo run --release --bin all-pats -- 250 >sols.json

# Compute the supercell matrices and a bunch of other info about the solutions.
# ("frac-lattice" in positions.json is the inverse of the supercell matrix)
cd ../hex-sols
python3 positions.py <../rust-sols/sols.json >positions.json

# Generate layers.yaml from them.
python3 make-inputs.py --params params/graphene.yaml <positions.json

# create MoS2 cells from position solutions, change "data" to position solution directory.
python3 MoS2.py data
'''
