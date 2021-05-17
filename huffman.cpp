#include <fstream>
#include <string>
#include <cstring>
#include <cstdint>
#include <string.h>
#include <map>
#include <iostream>
#include <iomanip>
#include <utility>
#include <vector>
#include <algorithm>

using namespace std;

void binDump(unsigned char i){
  for(int j = 0; j < 8; j++) std::cout << (i & (0x01 << (7-j) ) ? "1" : "0");
  std::cout << std::endl;
}

bool comparePair(const std::pair<char,int>& o1, const std::pair<char,int>& o2){
  return (o1.second > o2.second);
}


class BitSymbol{
  private:
    unsigned char length = 0;
    uint64_t data = 0;
    unsigned char iteratorPtr = 0;
  public:
    void add(bool bit){
      data |= (bit ? (uint64_t(1) << (63-length) ) : 0);
      length++;
    }
    void resetIterator(){
      iteratorPtr = 0;
    }
    bool hasNext(){
      return iteratorPtr < length;
    }
    bool getNext(){
      return data & (uint64_t(1) << (63 - (iteratorPtr++) ));
    }
    std::string getAsString() const {
      std::string ret;
      for(unsigned char i=0; i<length; i++){
        ret += (data & (uint64_t(1) << (63 - i)) ? "1" : "0");
      }
      return ret;
    }
    unsigned char getLength() const {
      return length;
    }
    bool equalsN(const BitSymbol &other, unsigned char n){
      uint64_t mask = 0;
      for(unsigned char i=0; i<n; i++){
        mask |= (uint64_t(1) << (63 - i));
      }
      return (this->data & mask) == (other.data & mask);
    }
    bool operator==(const BitSymbol o2) const {
      return (this->data == o2.data && this->length == o2.length);
    }
    bool operator<(const BitSymbol o2) const {
      return this->data < o2.data;
    }
    bool getLSB(){
      if(length==0) return false;
      return data & (uint64_t(1) << (64-length) );
    }
};

class BitStream{
  private:
    std::vector<unsigned char> bytes;
    unsigned char tmp = 0;
    char ptr = 0;
  public:
    static BitStream createFromString(const std::string& str, unsigned int startIndex=0){
      BitStream bs;
      for(unsigned int i=startIndex; i<str.length(); i++){
        bs.bytes.push_back(str[i]);
      }
      return bs;
    }
    void add(bool symbol){
      tmp |= (symbol ? (0x01 << (7-ptr) ) : 0);
      ptr++;
      if(ptr > 7){
        ptr = 0;
        bytes.push_back(tmp);
        tmp = 0;
      }
    }
    void addByte(unsigned char add){
      this->bytes.push_back(add);
    }
    void finalize(){
      while(ptr > 0) this->add(1);
    }
    std::vector<unsigned char> &getData(){
      return this->bytes;
    }
    std::string getAsString(){
      std::string ret;
      for(unsigned char b : this->bytes) ret += b;
      return ret;
    }
    unsigned int getLength() const {
      return bytes.size()*8 + ptr;
    }
    BitSymbol getSubBits(unsigned int start, unsigned int length){
      unsigned int startByte = start/8;
      unsigned int startBit = start%8;
      unsigned int endByte = (start+length-1)/8;
      unsigned int endBit = (start+length-1)%8;
      BitSymbol bs;
      for(unsigned int byteNum=startByte; byteNum<=endByte; byteNum++){
        if(byteNum >= this->bytes.size()) break;
        unsigned int bitStart = 0;
        if(byteNum == startByte) bitStart = startBit;
        unsigned int bitEnd = 7;
        if(byteNum == endByte) bitEnd = endBit;
        for(unsigned int bitNum=bitStart; bitNum<=bitEnd; bitNum++){
          bs.add(this->bytes[byteNum] & (0x01 << (7-bitNum)));
        }
      }
      return bs;
    }
    bool getSubBit(unsigned int index){
      unsigned int byteNum = index/8;
      unsigned int bitNum = index%8;
      if(byteNum >= this->bytes.size()) return false;
      return this->bytes[byteNum] & ( 0x01 << (7-bitNum) );
    }
};

class Huffman{
  private:
    static std::map<BitSymbol,char> generateSymbols(const std::vector<std::pair<char,int>> &sortedSymbolFrequencies){
      std::cout << "Count of different symbols: " << sortedSymbolFrequencies.size() << std::endl;
      std::vector<BitSymbol> bitSymbols;
      unsigned int charCount = 0;
      for(unsigned int i=0; i<sortedSymbolFrequencies.size(); i++){
        charCount += sortedSymbolFrequencies[i].second;
        BitSymbol bs;
        bitSymbols.push_back(bs);
      }

      int depth = 0;
      while(true && depth < 64){
        //int i=0; //only for DEBUG
        depth++;
        bool add = 0;
        unsigned int pi = 0;
        int compoundShare = 0;
        bool allFinal = true;
        bool lastLSB = 1;
        int bucketHalf;
        int bucketRelative;
        int bucketNumber = -1;
        for(const auto &p : sortedSymbolFrequencies){
          int share = p.second;
          compoundShare += share;
          //bool flipped = false; //only for DEBUG
          //bool finalSymbol = false; //only for DEBUG
          bool flagFirstInBucket = false;
          if(bitSymbols[pi].getLSB() != lastLSB){ //bucket changed
            bucketNumber++;
            lastLSB = bitSymbols[pi].getLSB();
            add = 0;
            flagFirstInBucket = true;
            //find next bucket
            int thisBucketCompound = share;
            for(unsigned int pj = pi+1; pj < sortedSymbolFrequencies.size(); pj++){
              if(bitSymbols[pj].getLSB() != bitSymbols[pi].getLSB()) break;
              thisBucketCompound += sortedSymbolFrequencies[pj].second;
            }
            bucketHalf = (thisBucketCompound)/2;
            bucketRelative = 0;
            if(share >= bucketHalf){
              bucketHalf = share;
            }
          }
          bucketRelative += share;
          if(bucketRelative >= bucketHalf && !flagFirstInBucket && !add){ //half bucket exceeded -> flip to 1
            add = 1;
            //flipped = true; //only for DEBUG
          }

          //check if it is done already
          bool differentThanNeighbors = true;
          if(pi > 0){ //has prev neighbor
            if(bitSymbols[pi].equalsN(bitSymbols[pi-1], bitSymbols[pi].getLength())) differentThanNeighbors = false;
          }
          if(pi < bitSymbols.size() - 1){ //has next neighbor
            if(bitSymbols[pi].equalsN(bitSymbols[pi+1], bitSymbols[pi].getLength())) differentThanNeighbors = false;
          }

          if(!differentThanNeighbors){ //not done - add a next bit
            bitSymbols[pi].add(add);
            allFinal = false;
          }/*else{ //this bit symbol is final
            finalSymbol = true; //only for DEBUG
          }*/
          /*std::cout << std::setw(2) << i++ << "| " << p.first << "(" << std::hex << ((int)p.first) << std::dec << ")"
                    << " -> " << std::setw(4) << p.second << "\t" << " (" << share << ") "
                    << compoundShare << "  B(" << bucketNumber << " : " << bucketRelative << " H" << bucketHalf << ")"
                    << "\t" << bitSymbols[pi].getAsString() << " " << (flipped ? "_" : "")
                    << (finalSymbol ? " F" : "") << std::endl; //only for DEBUG */
          pi++;
        }
        if(allFinal){
          //std::cout << " - - - ALL FINAL: WE ARE DONE - - - " << std::endl;
          break;
        }
      }

      std::map<BitSymbol,char> symbolSubstMap;
      for(unsigned int i = 0; i < sortedSymbolFrequencies.size(); i++){
        char c = sortedSymbolFrequencies[i].first;
        BitSymbol bitSymbol = bitSymbols[i];
        std::cout << std::setw(2) << std::dec << i << " | '" << c << "' (" << std::hex << ((int)c) << ") [" << std::dec << sortedSymbolFrequencies[i].second << "] -> " << bitSymbol.getAsString() << std::endl;
        symbolSubstMap[bitSymbol] = c;
      }

      return symbolSubstMap;
    }
  public:
    static std::pair<BitStream,std::map<BitSymbol,char>> strEncode(std::string in){
      unsigned int inLength = in.length();
      std::map<char,int> symbolFrequencies;
      unsigned int charCount = 0;
      for(std::string::iterator it=in.begin(); it!=in.end(); ++it){
        char c = *it;
        //std::cout << "in: '" << c << "' (" << std::hex << ((int)c) << ")" << std::endl;
        symbolFrequencies[c]++;
        charCount++;
      }
      std::vector<std::pair<char,int>> symbolsSort;
      //std::vector<BitSymbol> bitSymbols;
      for(const std::pair<char,int> &p : symbolFrequencies){
        symbolsSort.push_back(p);
        //BitSymbol bs;
        //bitSymbols.push_back(bs);
      }
      std::sort(symbolsSort.begin(), symbolsSort.end(), comparePair);

      //sorted, now generate symbols
      std::map<BitSymbol,char> symbolSubstMap = Huffman::generateSymbols(symbolsSort);
      int i = 0;
      for(const std::pair<BitSymbol,char> &p : symbolSubstMap){
        std::cout << std::setw(2) << i++ << " subst: " << p.first.getAsString() << " ~> '" << p.second << "'" << std::endl;
      }
      //we have all symbols generated

      BitStream bitStream;
      uint64_t progress = 0;
      std::cout << "Progress:" << std::endl;
      for(std::string::iterator it=in.begin(); it!=in.end(); ++it){
        //foreach char
        char c = *it;
        progress++;
        if(progress % 300 == 0){
          uint64_t percentage = uint64_t(10000) * progress / inLength;
          uint64_t percentageInt = percentage / 100;
          uint64_t percentageSub = percentage % 100;
          std::cout << "\r " << std::dec << percentageInt << "." << (percentageSub < 10 ? "0" : "") << percentageSub << "% ";
        }
        bool substFound = false;
        //find symbol
        for(const std::pair<BitSymbol,char> &symbolSubst : symbolSubstMap){
          if(symbolSubst.second != c) continue;
          //we have current index of current character
          //fetch our symbol:
          BitSymbol currentSymbol = symbolSubst.first;
          currentSymbol.resetIterator();
          while(currentSymbol.hasNext()){
            bitStream.add(currentSymbol.getNext());
          }
          substFound = true;
          break;
        }
        if(!substFound){ //error
          std::cout << "Error! Substitution for '" << c << "' (" << std::hex << ((int)c) << ") not found!" << std::endl;
        }
      }
      bitStream.finalize();


      std::vector<unsigned char> rawBytes = bitStream.getData();
      unsigned int outLength = rawBytes.size();
      /*for(unsigned char b : rawBytes){
        std::cout << "-- ";
        binDump(b);
      }*/
      std::cout << "Input length: " << inLength << ", output length: " << outLength << std::endl;
      return std::make_pair(bitStream, symbolSubstMap);
    }

    static std::string decode(BitStream enc, std::map<BitSymbol,char> symbolSubstMap){
      unsigned int stringPos = 0;
      unsigned int encLength = enc.getLength();
      unsigned int progress = 0;
      std::string dec = "";
      std::cout << "Progress:" << std::endl;
      while(true){
        progress++;
        if(progress % 300 == 0){
          uint64_t percentage = uint64_t(10000) * stringPos / encLength;
          uint64_t percentageInt = percentage / 100;
          uint64_t percentageSub = percentage % 100;
          std::cout << "\r " << std::dec << percentageInt << "." << (percentageSub < 10 ? "0" : "") << percentageSub << "% ";
        }
        bool foundSymbol = false;
        for(const std::pair<BitSymbol,char> &p : symbolSubstMap){
          unsigned int symbolLen = (unsigned int) p.first.getLength();
          BitSymbol symbol = enc.getSubBits(stringPos, symbolLen);
          if(p.first == symbol){ //symbol found
            dec += p.second;
            stringPos += symbolLen;
            //std::cout << "Matched char '" << p.second << "'!" << std::endl;
            foundSymbol = true;
            break;
          }
        }
        if(!foundSymbol){
          std::cout << "Symbol not matched! stringPos=" << stringPos << std::endl;
          break;
        }
        if(enc.getLength() <= stringPos){
          std::cout << "End reached! enc.getLength()=" << enc.getLength() << " <= stringPos=" << stringPos << std::endl;
          break;
        }
      }
      return dec;
    }

    static std::string serialize(BitStream bitStream, std::map<BitSymbol,char> symbolSubstMap){
      std::string ret = "";
      ret += ((unsigned char)0xAD); //header part 1
      ret += ((unsigned char)0xBD); //header part 2
      ret += ((unsigned char)0x01); //version 1
      ret += ((unsigned char)(symbolSubstMap.size())); //number of substitution symbol
      //table of subst symbols
      BitStream symbolMapStream;
      for(std::pair<BitSymbol,char> p : symbolSubstMap){
        //char (8b), symbolSize (6b), symbol (Xb)
        for(int i=0; i<8; i++)
          symbolMapStream.add(p.second & (0x01 << (7-i) ));
        unsigned char symbolSize = p.first.getLength();
        for(int i=2; i<8; i++)
          symbolMapStream.add(symbolSize & (0x01 << (7-i) ));
        p.first.resetIterator();
        while(p.first.hasNext())
          symbolMapStream.add(p.first.getNext());
      }
      symbolMapStream.finalize();
      const std::vector<unsigned char> &symbols =  symbolMapStream.getData();
      for(const unsigned char &c : symbols){
        ret += c;
      }
      //raw compressed data
      //maybe 8 B of compressed data length
      bitStream.finalize(); //just to be sure
      const std::vector<unsigned char> &symbolCharacters =  bitStream.getData();
      for(const unsigned char &c : symbolCharacters){
        ret += c;
      }
      return ret;
    }

    static std::pair<BitStream,std::map<BitSymbol,char>> deserialize(const std::string& serial){
      if(serial.length() < 4) return std::make_pair<BitStream,std::map<BitSymbol,char>>(BitStream(),std::map<BitSymbol,char>()); //invalid header
      if(((unsigned char)serial[0]) != 0xAD || ((unsigned char)serial[1]) != 0xBD) return std::make_pair<BitStream,std::map<BitSymbol,char>>(BitStream(),std::map<BitSymbol,char>()); //invalid format
      if(serial[2] != 0x01) return std::make_pair<BitStream,std::map<BitSymbol,char>>(BitStream(),std::map<BitSymbol,char>()); //not supported format version
      // /\ these will be exceptions
      unsigned char symbolSubstMapSize = (unsigned char)serial[3];
      std::map<BitSymbol,char> symbolSubstMap;
      BitStream bitStream = BitStream::createFromString(serial, 4);
      unsigned int bitCount = 0;
      //parse symbol substitution map
      for(unsigned char i=0; i<symbolSubstMapSize; i++){
        //char (8b), symbolSize (6b), symbol (Xb)
        char c = 0;
        for(int j=0; j<8; j++){
          c |= (bitStream.getSubBit(bitCount++) ? (0x01 << (7-j) ) : 0);
        }
        unsigned char symbolSize = 0;
        for(int j=0; j<6; j++){
          bool got = bitStream.getSubBit(bitCount++);
          symbolSize |= (got ? (0x01 << (5-j) ) : 0);
        }
        if(symbolSize == 0) symbolSize = 64;
        BitSymbol bitSymbol;
        for(unsigned char j=0; j<symbolSize; j++){
          bitSymbol.add(bitStream.getSubBit(bitCount++));
        }
        symbolSubstMap[bitSymbol] = c;
        std::cout << "D" << std::dec << std::setw(2) << (unsigned int)i << " subst: " << bitSymbol.getAsString() << " ~> '" << c << "'" << std::endl;
      }
      //discard padding to next byte
      while(bitCount % 8 != 0) bitCount++;
      //copy data
      BitStream dataBitStream;
      unsigned int dataStart = 4 + bitCount/8;
      for(unsigned int i=dataStart; i<serial.length(); i++){
        dataBitStream.addByte(serial[i]);
      }
      return std::make_pair(dataBitStream, symbolSubstMap);
    }
};


class File{
  private:
    const char* filename;
  public:
    File(const char* f){
      this->filename = f;
    }
    string read(){
      string c = "";
      string s;
      fstream fs;
      fs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
      try{
        fs.open(this->filename, ios_base::in);
        while(getline(fs, s)){
          c += s+"\n";
        }
        if(c.length()){
          c.pop_back();
        }
        fs.close();
      }catch(const std::fstream::failure &e){
        std::cerr<<"Error in read: "<<e.what()<<std::endl;
        cerr << "Error: " << strerror(errno) << std::endl;
      }
      return c;
    }
    bool write(const string& w){
      fstream fs;
      fs.exceptions(std::ifstream::failbit | std::ifstream::badbit);
      try{
        fs.open(this->filename , ios::trunc|ios::out);
        fs<<w;
        fs.close();
      }catch(const std::fstream::failure &e){
        std::cerr<<"Error in write: "<<e.what()<<std::endl;
        cerr << "Error: " << strerror(errno) << std::endl;
        return false;
      }
      return true;
    }
    void append(const string& w) {
      fstream fs;
      fs.open(this->filename, ios::out | ios::app);
      fs<<w;
      fs.close();
    }
};

int main(int argc, char** argv){
    if(argc != 3){
        std::cout << "Usage: " << argv[0] << " file_in file_out" << std::endl;
        return 0;
    }
    std::string filenameIn = std::string(argv[1]);
    std::string filenameOut = std::string(argv[2]);
    std::cout << "Encoding file " << filenameIn << " outputing to file " << filenameOut << std::endl;

    //read file
    File fileIn(filenameIn.c_str());
    std::string in = fileIn.read();
    //encrypt - different section
    std::pair<BitStream,std::map<BitSymbol,char>> out = Huffman::strEncode(in);

    //TESTING
    std::string serial = Huffman::serialize(out.first, out.second);
    std::cout << "serial length=" << serial.length() << " B" << std::endl;
    std::pair<BitStream,std::map<BitSymbol,char>> unserial = Huffman::deserialize(serial);
    std::cout << "unserial size=" << unserial.first.getLength() << " b; " << unserial.second.size() << " size" << std::endl;

    std::cout << "Trying to decode!" << std::endl;
    std::string dec = Huffman::decode(out.first, out.second);
    std::string dec2 = Huffman::decode(unserial.first, unserial.second);
    //std::cout << dec;

    //writing output file
    std::cout << "1f: " << filenameOut << std::endl;
    File fileOut(filenameOut.c_str());
    fileOut.write(out.first.getAsString());

    std::string filenameOutDec = filenameOut+".dec";
    std::string filenameOutDec2 = filenameOut+".dec2";
    std::cout << "2f: " << filenameOutDec << " = " << filenameOutDec.c_str() << std::endl;
    std::cout << "3f: " << filenameOutDec2 << " = " << filenameOutDec2.c_str() << std::endl;
    File fileOut2(filenameOutDec.c_str());
    fileOut2.write(dec);
    File fileOut3(filenameOutDec2.c_str());
    fileOut3.write(dec2);
    return 0;
}
