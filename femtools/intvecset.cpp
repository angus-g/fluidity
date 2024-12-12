#include <vector>
#include <set>
#include <list>

typedef std::set<std::vector<int> > intvecset;
typedef std::list<intvecset> intveclist;
intveclist IVS;

extern "C"
{
  void intvec_create_set(int* idx);
  void intvec_is_present(int* idx, int* arr, int* size, int* success);
  void intvec_clear_set(int* idx);
  void intvec_destroy_set(unsigned int* idx);
}

void intvec_create_set(int* idx)
{
  intvecset S;
  IVS.push_back(S);
  *idx = IVS.size();
}

void intvec_is_present(int* idx, int* arr, int* size, int* success)
{
  std::vector<int> v(*size);
  std::pair<intvecset::iterator,bool> stat;
  int i;
  intvecset S;

  intveclist::iterator j;
  int k;
  j = IVS.begin();
  for (k = *idx; k > 0; k--)
    j++;
  S = *j;

  for (i = 0; i < *size; i++)
    v[i] = arr[i];

  stat = S.insert(v);
  *success = stat.second;
}

void intvec_clear_set(int* idx)
{
  intvecset S;
  intveclist::iterator j;
  int k;
  j = IVS.begin();
  for (k = *idx; k > 0; k--)
    j++;
  S = *j;
  S.clear();
}

void intvec_destroy_set(unsigned int* idx)
{
  intvecset S;
  intveclist::iterator j;
  int k;
  j = IVS.begin();
  for (k = *idx; k > 0; k--)
    j++;
  S = *j;
  S.clear();
  if (IVS.size() == *idx)
  {
    IVS.pop_back();
  }
}
