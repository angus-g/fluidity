#include <set>

#ifdef DOUBLEP
#define REAL double
#else
#define REAL float
#endif

typedef std::set<int> eleset;
eleset E;

extern "C"
{
  void ele_add_to_set(int element);
  void ele_get_size(int* size);
  void ele_fetch_list(int* arr);
  void ele_get_ele(int i, int* ele);
}

void ele_add_to_set(int element)
{
  E.insert(element);
}

void ele_get_size(int* size)
{
  *size = E.size();
}

void ele_fetch_list(int* arr)
{
  int pos;

  pos = 0;
  for (eleset::const_iterator i = E.begin(); i != E.end(); i++)
  {
    arr[pos++] = *i;
  }

  E.clear();
}

void ele_get_ele(int i, int* ele)
{
  int pos;

  pos = 0;
  for (eleset::const_iterator j = E.begin(); j != E.end(); j++)
  {
    pos = pos + 1;
    if (pos == i)
    {
      *ele = *j;
      return;
    }
  }
}
