// Scheduling Arduino Code (run every x seconds)
// http://hwhacks.com/2016/05/08/scheduling-arduino-code-run-every-x-seconds/#:~:text=Arduino%20code%20by%20nature%20runs,halt%20the%20program%20before%20continuing.

byte led = 2;

typedef struct t  {
    unsigned long tStart;
    unsigned long tTimeout;
};

//Tasks and their Schedules.
t t_func1 = {0, 1}; //Need one function to run every 10 microseconds.  
//t t_func2 = {0, 2000}; //Run every 2 seconds. 

bool tCheck (struct t *t ) {
  if (millis() > t->tStart + t->tTimeout) return true;    
}

void tRun (struct t *t) {
//change this to micros()    
    t->tStart = micros();
}

void setup (void) {
  //Arduino setup.
  pinMode(led, OUTPUT);
}

void loop (void) {
    if (tCheck(&t_func1)) {
      func1();
      tRun(&t_func1);
    }
    
//    if (tCheck(&t_func2)) {
//      func2();
//      tRun(&t_func2);
 //   }
}

void func1 (void) {
  //This executes every 6s.
  digitalWrite(led, HIGH); 
  digitalWrite(led, LOW);
}

//void func2 (void) {
  //This executes every 2 seconds.
//}
