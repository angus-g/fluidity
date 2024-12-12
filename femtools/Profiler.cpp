/* Copyright (C) 2010- Imperial College London and others.

   Please see the AUTHORS file in the main source directory for a full
   list of copyright holders.

   Dr Gerard J Gorman
   Applied Modelling and Computation Group
   Department of Earth Science and Engineering
   Imperial College London

   g.gorman@imperial.ac.uk

   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with this library; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
   USA
*/

#include "Profiler.h"

using namespace std;

Profiler::Profiler(){
}

Profiler::~Profiler(){}

double Profiler::get(const std::string &key) const{
  double time = timings.find(key)->second.second;
#ifdef HAVE_MPI
  int init_flag;
  MPI_Initialized(&init_flag);
  if(init_flag){
    double gtime;
    MPI_Reduce(&time, &gtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    return gtime;
  }
#endif

  return time;
}

void Profiler::print() const{
  bool print = true;
  double val;
#ifdef HAVE_MPI
  int MyRank, init_flag;
  MPI_Comm_rank(MPI_COMM_WORLD, &MyRank);
  MPI_Initialized(&init_flag);
  print = !init_flag || MyRank == 0;
#endif
  for(map< string, pair<double, double> >::const_iterator it=timings.begin();it!=timings.end();++it){
    val = get(it->first);
    if ( print ) {
      cout<<it->first<<" :: "<<val<<endl;
    }
  }
}

void Profiler::tic(const std::string &key){
  timings[key].first = wall_time();
}

void Profiler::toc(const std::string &key){
  timings[key].second += wall_time() - timings[key].first;
}

double Profiler::wall_time() const{
#ifdef HAVE_MPI
  return MPI_Wtime();
#else
  return 0.0;
#endif
}

void Profiler::zero(){
  for(map< string, pair<double, double> >::iterator it=timings.begin();it!=timings.end();++it){
    it->second.second = 0.0;
  }
}

void Profiler::zero(const std::string &key){
  timings[key].second = 0.0;
}

int Profiler::minorpagefaults(){
  int faults = -99;
#ifdef HAVE_LIBNUMA
  getrusage(RUSAGE_SELF, &flprofiler.usage);
  faults = flprofiler.usage.ru_minflt;
#endif
  return faults;
}

int Profiler::majorpagefaults(){
  int faults = -99;
#ifdef HAVE_LIBNUMA
  getrusage(RUSAGE_SELF, &flprofiler.usage);
  faults = flprofiler.usage.ru_majflt;
#endif
  return faults;
}


int Profiler::getresidence(void *ptr){
  int residence=-99;
#ifdef HAVE_LIBNUMA
  int mode;
  size_t page_size = getpagesize();
  size_t page_id = (size_t)ptr/page_size;
  /* round memory address down to start of page */
  void *start_of_page =  (void *)(page_id*page_size);

  /* If flags  specifies  both MPOL_F_NODE and MPOL_F_ADDR,
   * get_mempolicy() will return the node ID of the node on
   * which the address of the start of the page is allocated
   * into the location pointed to by mode
   */
  unsigned long flags = MPOL_F_NODE|MPOL_F_ADDR;

  // if(get_mempolicy(&mode, NULL, 0, start_of_page, flags)){
  //   perror("get_mempolicy()");
  // }
  // residence = mode;
#endif
  return residence;
}

// Opaque instances of profiler.
Profiler flprofiler;

// Fortran interface
extern "C" {
  void cprofiler_get(const char *key, const int key_len, double *time){
    *time = flprofiler.get(string(key, key_len));
  }

  void cprofiler_tic(const char *key, const int key_len){
    flprofiler.tic(string(key, key_len));
  }

  void cprofiler_toc(const char *key, const int key_len){
    flprofiler.toc(string(key, key_len));
  }

  void cprofiler_zero(){
    flprofiler.zero();
  }

  void cprofiler_minorpagefaults(int *faults){
    *faults = flprofiler.minorpagefaults();
  }

  void cprofiler_majorpagefaults(int *faults){
    *faults = flprofiler.majorpagefaults();
  }

  void cprofiler_getresidence(void *ptr, int *residence){
    *residence = flprofiler.getresidence(ptr);
  }
}
