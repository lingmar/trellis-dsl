# DNA sequencing

## Basic definitions
* A *DNA string* (or *genome*): a sequence of nucleotides. 
* A *nucleotide* (or *base*): one of four molecule in a DNA string, abbreviated
  A, C, G, and T. Each nucleotide is mapped bi-directionally to integers via as
  A=0, C=1, G=2, and T=3. Let `base2idx` be the function that maps bases to
  integers and `idx2base` be its inverse.
  
* A *k-mer*: a DNA string of length `k`. Typically, `k` is short (3-7), and the
  k-mer is a substring of a longer DNA string. A set of k-mers are
  bi-directionally mapped to integers by interpreting them as base 4 numbers.
  For example, the 3-mer
  ```
  ACT
  ```
  is mapped into
  ```
  0 * 4^0 + 1 * 4^1 + 3 * 4^2 = 49
  ```
  
  Let `last` return the last (rightmost) base of a k-mer. For example,
  `last(ACT) = T`.
  
## Explicit Duration Model

### Hidden states
Each hidden state is represented by a pair (`kmer`, `depth`). The `kmer`
represents the k-mer in the genome currently being measured. The `depth`
represents the number of time steps the `kmer` will stay in the machine.

The total number of states is `4^k*D` where `k` is the length of the k-mers and
`D` is the maximum depth.

A state transition from a state (`kmer1`, `depth1`) to another state (`kmer2`,
`depth2`) is in one of these cases:
1. `depth1=0` (assuming indexing from `0`). Then the genome shifts one
  base. That is, `kmer2` shifts one base relative to `kmer1`. For example, the
  3-mer `ATC` may shift into `TCA`, `TCC`, `TCG`, or `TCT`.
2. `depth1` is smaller than the maximum depth but larger than `0`. Then `kmer2`
  is equal to `kmer1` and `depth2` is one smaller than `depth1`
3. `depth1` is equal to the maximum depth. Then `depth2` *may* decrease one step
  but does not have to. The k-mer does not shift, and `kmer2` is thus equal
  to `kmer1`.

### Raw input
  The raw inputs are given via files in fast 5 formats.
  * A model file, for example `synthetic_dna/model/Params3_true.h5` contains model
    parameters as well as data that has been found via training. It has the
    following layout:
  
  ```
  .
  ├── Parameters
  │   └── Normalization
  │       ├── HighClip     # ?
  │       ├── LowClip      # ?
  │       └── SignalLevels # Number of signal levels
  └── Tables
      ├── DurationProbabilities    # Array of floats of length D
      ├── ObservationProbabilities # Matrix of floats of size 4^k*SignalLevels
      │   └── TailFactor           # Int
      └── TransitionProbabilities  # Matrix of floats of size 4^k*4
  ```
  The meaning of the probability tables are:
  * `ObservationProbabilities[i,s] = P(s | i)`, i.e. the probability of
    observing signal `s` given the `i`:th k-mer (regardless of depth).
  * `DurationProbabilities[d]`: the probability of the current k-mer staying for
     `d` time steps.
  * `TransitionProbabilities[i,b]`: the probability of the `i`:th k-mer shifting
    so that its "newest" base becomes `b`.

* A signal file, for example `synthetic_dna/synthetic_dna/data/signals100_1.fast5`:
  ```
  ├── <key 1>          # Instance-dependent key
  │   └── Raw
  │       └── Signal   # Array of ints ranging from 0 to SignalLevels
  ...
  └── <key n>
      └── Raw
          └── Signal
  ```
  
* A mapped file, for example `synthetic_dna/data/mapped_reads100_1.hdf5`:
  ```
  ├── <key 1>            # Same keys as in the signal file
  │   └── Reads
  │       └── Reference  # Array of ints representing bases (0=A,1=C,2=G,3=T)
  └── key2
      └── Reads
          └── Reference
  ```

### Model parameters
* `k`: the number of bases in a k-mer in the hidden states. Found by inspecting
  the size of the `ObservationProbabilities` table.
* `D`: the maximum depth in a hidden state. Found by inspecting the size of the
  `DurationProbability` vector.
* `SignalLevels`: the maximum signal value.

### Input signal
The input signal is a vector of integers ranging from `0` to `SignalLevels`,
given by a `Signal` entry from a signal file.

### Observation probabilities
The probability of observing a signal `s` given a hidden state (`kmer`, `depth`)
is `ObservationProbabilities[base2idx(kmer), s]`.

### Transition probabilities
The probability of transitioning from (`kmer1`, `depth1`) into (`kmer2`,
`depth2`) follow the three cases for transitions listed above:
1. `depth1 = 0`: The probability of shifting from `kmer1` to `kmer2` multiplied
   by the probability of `kmer2` staying for `depth2` time steps:
   
   ```
   TransitionProbabilities[base2idx(kmer1), last(kmer2)] * DurationProbabilities[depth2]`
   ```
   
2. `0 <= depth1 <= D`: The only state that (`kmer1`, `depth1`) can transition
   into is (`kmer1`, `depth1` - 1). For that state the transition probability is
   1, for all others it is 0.
   
3. `depth = D`: The state (`kmer1`, `depth1`) transitions into itself with
   probability `TailFactor`, and into (`kmer1`, `depth1` - 1) with probability
   `1-TailFactor`.

## Explicit Duration Model
TODO
