#include <stdlib.h>
// #include "test.h"
struct TestStructure {/** */
  int length;
  char title;
};

int test(int ntest,struct TestStructure **test_ptr_ptr){
  for (int i=0;i<ntest;i++){
      printf("%d\n",i);
      printf("%d\n", &test_ptr_ptr[0][i]);
      printf("%d\n", test_ptr_ptr[0][i]);
      if (test_ptr_ptr[i] == 0){
        printf("foo\n");
      }
      printf("length - %d\n", test_ptr_ptr[0][i].length);
      if (test_ptr_ptr[0][i].length == 0){
        printf("bar\n");
      }
  }
  return 1;
}

struct TestStructure *outputTestStructure(){
  struct TestStructure *p1 = malloc(10 * sizeof(struct TestStructure *));
  struct TestStructure *a1 = p1;
  for(int i=0;i<10;i++){
      a1[i] = (struct TestStructure){i,'b'};
  }
  return a1;
}

main(){
  struct TestStructure *a1 = outputTestStructure();
  struct TestStructure **test_ptr_ptr = &a1;
  test(10,test_ptr_ptr);
}