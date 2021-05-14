# cpp-huffman-coding-compression
Simple lossless compression algorithm for arbitrary data

## Huffman coding
Huffman coding is a technique where you sort symbols in a file by frequency and you assign the more frequent symbols the shorter symbols and the less frequent symbols the longer symbols. An average text for example contains letters 'e' and 'a' a lot more more than letters 'x' and 'w'. So you can transform all 'e' letters into 3-bit sequence and all 'x' letters into for example 20-bit sequence. Because 'e' is a lot more frequent than 'x', the resulting bit stream will probably be a lot smaller than the original text without any loss of information. When decoding, you just substitute back the right letters in place of their bit symbols.
