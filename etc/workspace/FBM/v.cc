#include <boost/version.hpp>
#include <boost/regex.hpp>
#include <iostream>

using namespace std;

#define EACHA(a) for (size_t i=0; i<sizeof(a)/sizeof((a)[0]); i++)
#define EACHV(v) for (size_t i=0; i<(v).size(); i++)
int main()
{
  vector<string> v;
  cout << BOOST_VERSION << endl;
  cout << BOOST_VERSION/100000 << endl;
  cout << BOOST_VERSION/100%1000 << endl;
  cout << BOOST_VERSION%100 << endl;

  v.push_back("hoge");
  v.push_back("fuga");
  v.push_back("piyo");
  v.push_back("foo");
  v.push_back("bar");
  v.push_back("baz");
  EACHA(v){
    cout << "each a\t" << v[i] << endl;
  }
  EACHV(v){
    cout << "each v\t" << v[i] << endl;
  }
}
