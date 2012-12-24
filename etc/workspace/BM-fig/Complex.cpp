#include <cppunit/TestCase.h>
#include <cppunit/TestResult.h>
#include <cppunit/TestResultCollector.h>
#include <cppunit/TestCaller.h>
#include <cppunit/CompilerOutputter.h>
#include <stdexcept>

class Complex
{
  friend bool operator==(const Complex& a, const Complex& b);
  
  double real, imaginary;
public:
  Complex(double r, double i=0)
    : real(r), imaginary(i){}

  Complex& operator+(const Complex& rhs)
  {
    Complex lhs(0.0);
    lhs.real = this->real + rhs.real;
    lhs.imaginary = this->imaginary + rhs.imaginary;
    return lhs;
  }
};

bool operator==(const Complex &a, const Complex &b)
{
  return a.real == b.real && a.imaginary == b.imaginary;
}


//class ComplexNumberTest : public CppUnit::TestCase
//{
//public:
//  ComplexNumberTest( std::string name ) : CppUnit::TestCase( name ){}
//
//  void runTest(){
//    CPPUNIT_ASSERT( Complex(10,1) == Complex(10,1));
//    CPPUNIT_ASSERT(!(Complex(1,1) == Complex(2,2)));
//  }
//};

class ComplexNumberTest : public CppUnit::TestFixture
{
private:
  Complex *m_10_1, *m_1_1, *m_11_2;
public:
  void setUp()
  {
    m_10_1 = new Complex( 10, 1 );
    m_1_1  = new Complex(  1, 1 );
    m_11_2 = new Complex( 11, 2 );
  }

  void tearDown()
  {
    delete m_10_1;
    delete m_1_1;
    delete m_11_2;
  }


  void testEquality()
  {
    CPPUNIT_ASSERT( *m_10_1 == *m_10_1  );

    CPPUNIT_ASSERT( *m_10_1 == *m_1_1  );
    CPPUNIT_ASSERT(!(*m_10_1 == *m_11_2));
    CPPUNIT_ASSERT( *m_11_2 == *m_1_1  );
  }

  void testAddition()
  {
    CPPUNIT_ASSERT( *m_10_1 + *m_1_1 == *m_11_2 );
  }

};


//class Complex
int main(int argc, char** argv)
{


  CppUnit::TestCaller<ComplexNumberTest> 
    test("testEquality", &ComplexNumberTest::testEquality );

  CppUnit::TestResult controller;
  CppUnit::TestResultCollector result;
  controller.addListener( &result );
  try{
    test.run( &controller );

    // Print test in a compiler comatible format.
    CppUnit::CompilerOutputter outputter( &result, std::cerr );
    outputter.write();
  }
  catch( std::invalid_argument &e ){
    std::cerr << std::endl
	      << "Error: " << e.what()
	      << std::endl;
  }
  return result.wasSuccessful() ? 0 : 1;
}
  

