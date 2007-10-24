
#include <stdlib.h>
#include <stdio.h>
#include <dlfcn.h>

/* Note that we don't have to include "libhello.h".
   However, we do need to specify something related;
   we need to specify a type that will hold the value
   we're going to get from dlsym(). */

/* The type "simple_demo_function" describes a function that
   takes no arguments, and returns no value: */

typedef void (*simple_demo_function)(void);
typedef void (*simple_open)(char *);


int main(void) {
 const char *error;
 void *module;
 simple_demo_function demo_function;
 simple_open open_function;
 /* Load dynamically loaded library */
 module = dlopen("perform.so", RTLD_LAZY);
 if (!module) {
   fprintf(stderr, "Couldn't open perform.so: %s\n",
           dlerror());
   exit(1);
 }

 
 open_function = dlsym(module, "readFile");
 if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find readFile: %s\n", error);
   exit(1);
 }

 /* Now call the function in the DL library */
 (*open_function)("/scratch/L.cir");
 demo_function = dlsym(module, "MEprint");
 if ((error = dlerror())) {
   fprintf(stderr, "Couldn't find MEprint: %s\n", error);
   exit(1);
 }
 (*demo_function)();

 /* All done, close things cleanly */
 dlclose(module);
 return 0;
}
