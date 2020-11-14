typedef struct t  {
    unsigned long tStart;
    unsigned long tTimeout;
};

typedef struct square_wave {
    unsigned long period;
    unsigned long phase;
    unsigned long duty;
    byte polarity;
    int counter; //this is what Dr. Dahl called 'clock'
    int pin;
    int state;
};


//pin = 37 
struct square_wave wave1 = {10,2,50,HIGH,0,37,0}; 


//Tasks and their Schedules.
t t_func1 = {0, 10}; 

bool tCheck (struct t *t ) {
  if (micros() > t->tStart + t->tTimeout) return true;    
}

void tRun (struct t *t) {
  t->tStart = micros();
}

void update_wave(square_wave wave1) {
if (wave1.state==0 and wave1.counter > wave1.phase){
  wave1.state = 1;
}
else if (wave1.state==1 and wave1.counter > (wave1.phase + wave1.duty)){
  wave1.state = 2;
}
else if (wave1.state==2 and wave1.counter > wave1.period){
     wave1.state = 0;
     wave1.counter = 0;
}

wave1.counter++;

if (wave1.state==1) {
  digitalWrite(wave1.pin,wave1.polarity); 
}
else { 
  digitalWrite(wave1.pin, not wave1.polarity); 
}

}

void setup (void) {
  //Arduino setup.
  pinMode(wave1.pin, OUTPUT);
 
}

void loop (void) {
  if (tCheck(&t_func1)) {
    func1();
    tRun(&t_func1);
    }
}

// the square wave lives in here, and eventually so too will the trigger fan in / fan out. 
void func1 (void) {
  update_wave(wave1);
}
  
