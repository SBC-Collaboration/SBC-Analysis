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

unsigned long DELAY_TIME = 10; // 10 us
unsigned long delayStart = 0; // the time the delay started
//bool delayRunning = false; // true if still waiting for delay to finish

//bool ledOn = false; // keep track of the led state


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
  digitalWrite(wave1.pin, LOW); // turn led off
//  ledOn = false;

  // start delay
  delayStart = micros();
//  delayRunning = true;
} 

void loop (void) {
  while ((micros() - delayStart) <= DELAY_TIME) {
    delayStart += DELAY_TIME; // this prevents drift in the delays
    update_wave(wave1);
  }
}


  
