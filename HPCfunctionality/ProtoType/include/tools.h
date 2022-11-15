#include <Eigen/Dense>


std::vector<double> collapseVectorStd(std::vector<double> const & vin, sbn::SBNconfig const & conf){
   // All we want is a representation with the subchannels added together
   std::vector<double> cvec(conf.num_bins_total_compressed, 0.0);
   for (int d=0; d<conf.num_detectors;++d) {
      size_t offset_in(0), offset_out(0);
      for (int i=0; i<conf.num_channels; i++) {
          size_t nbins_chan = conf.num_bins[i];
          for (int j=0; j<conf.num_subchannels[i]; j++) {
             size_t first_in   = d*conf.num_bins_detector_block            + offset_in;
             size_t first_out  = d*conf.num_bins_detector_block_compressed + offset_out;
             std::transform (
                   cvec.begin() + first_out, cvec.begin() + first_out + nbins_chan, 
                   vin.begin()  + first_in,  cvec.begin() + first_out,
                   std::plus<double>());
             offset_in +=nbins_chan;
          }
          offset_out += nbins_chan;
      }
   }
   return cvec;
}

// TODO can this be simplified?
std::vector<double> collapseVector(std::vector<double> const & vin, sbn::SBNconfig const & conf){
   std::vector<double> cvec;
   cvec.reserve(conf.num_bins_total_compressed);
   //std::vector<double> cvec2(conf.num_bins_total_compressed, 0.0);

   for(int im = 0; im < conf.num_modes; im++){
       for(int id =0; id < conf.num_detectors; id++){
           int edge = id*conf.num_bins_detector_block + conf.num_bins_mode_block*im;
           for(int ic = 0; ic < conf.num_channels; ic++){
               int corner=edge;
               for(int j=0; j< conf.num_bins.at(ic); j++){
                   double tempval=0;
                   for(int sc = 0; sc < conf.num_subchannels.at(ic); sc++){
                        tempval += vin.at(j+sc*conf.num_bins.at(ic)+corner);
                        edge +=1;
                   }
                   cvec.push_back(tempval);
               }
           }
       }
   }
   return cvec;
}

// Split input vector into pieces of similar size
template <typename T>
std::vector<std::vector<T>> splitVector(const std::vector<T>& vec, size_t n) {
    std::vector<std::vector<T>> outVec;
    size_t length = vec.size() / n;
    size_t remain = vec.size() % n;
    size_t begin = 0;
    size_t end = 0;

    for (size_t i = 0; i < std::min(n, vec.size()); ++i) {
        end += (remain > 0) ? (length + !!(remain--)) : length;
        outVec.push_back(std::vector<T>(vec.begin() + begin, vec.begin() + end));
        begin = end;
    }
    return outVec;
}

// Concatenate vector of vectors as vector
std::vector<double> asVector(std::vector<std::vector<double> > v_in) {
    std::vector<double> allSpectra;
    allSpectra.reserve(v_in.size()*v_in[0].size());
    for (auto temp : v_in) allSpectra.insert(allSpectra.end(), temp.begin(), temp.end());
    return allSpectra;
}

std::vector<double> asVector(TMatrixT<double> const & M) {
    const double *pData = M.GetMatrixArray();
    std::vector<double> vData;
    vData.assign(pData, pData + M.GetNoElements());
    return vData;
}

template<typename T>
std::vector<T> myslice(std::vector<T> const &v, int m, int n)
{
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;
    std::vector<T> vec(first, last);
    return vec;
}

std::vector< std::vector<double> > sliceVector(std::vector<double> const & input, int nPoints) {
   std::vector< std::vector<double> > test;
   test.reserve(nPoints);
   int nBins = input.size()/nPoints;

   std::vector<double> work;
   work.reserve(nBins);

   for (int i=0;i<nPoints;++i) {
      work=myslice(input, i*nBins, (i+1)*nBins-1);
      test.push_back(work);
   }
   return test;
}

// Properly release memory
void releaseVec(std::vector<double> & vec) {
    vec.clear();
    vec.shrink_to_fit();
}
