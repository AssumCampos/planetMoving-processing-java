// Authors:  //<>//
// Maria Assumpció Campos Martínez 

// Problem description:
// A planet moving around a star

// Differential equations:
// F = f1
// f1 = G*M*m/r^2 

// Definitions:

enum IntegratorType 
{
  NONE,
  EXPLICIT_EULER, 
  SIMPLECTIC_EULER, 
  HEUN, 
  RK2, 
  RK4
}

// Parameters of the numerical integration:

final boolean REAL_TIME = true;
float SIM_STEP = 20000;   // Simulation time-step (s)
final float TIME_AGES = 1E6;
IntegratorType integrator = IntegratorType.EXPLICIT_EULER;   // ODE integration method
String integrador = "Euler";
String paso_sim = str(SIM_STEP);

// Display values:

final boolean FULL_SCREEN = true;
final int DRAW_FREQ = 50;   // Draw frequency (Hz or Frame-per-second)
//int DISPLAY_SIZE_X = 1000;   // Display width (pixels)
//int DISPLAY_SIZE_Y = 1000;   // Display height (pixels)

int DISPLAY_SIZE_X = 1920;   // Display width (pixels)
int DISPLAY_SIZE_Y = 1080;   // Display height (pixels)

// Draw values:

final int [] BACKGROUND_COLOR = {200, 200, 255};
final int [] REFERENCE_COLOR = {0, 255, 0};
final int [] OBJECTS_COLOR = {255, 0, 0};
final float OBJECTS_SIZE = 10E9;   // Size of the objects (m)
final float PIXELS_PER_METER = 2E-9;   // Display length that corresponds with 1 meter (pixels)
final PVector DISPLAY_CENTER = new PVector(0.0, 0.0);   // World position that corresponds with the center of the display (m)

// Parameters of the problem:

final float Me = 1.989E30;   // Star mass (kg)
final float Mp = 5.972E24; // Planet mass (kg)

final float Gc = 6.67430E-11;   // Gravity universalconstant (N * (m*m)/(kg*kg))
final PVector G = new PVector(0.0, -Gc);   // Acceleration due to gravity (m/(s*s))


// Time control:
int _lastTimeDraw = 0;   // Last measure of time in draw() function (ms)
float _deltaTimeDraw = 0.0;   // Time between draw() calls (s)
float _simTime = 0.0;   // Simulated time (s)
float _elapsedTime = 0.0;   // Elapsed (real) time (s)

// Output control:

PrintWriter _output;
final String FILE_NAME = "data.txt";

// Auxiliary variables:

float _energy;   // Total energy of the particle (J)

// Variables to be solved:

PVector _s = new PVector();   // Position of the planet (m)
PVector _v = new PVector();   // Velocity of the planet (m/s)
PVector _a = new PVector();   // Accleration of the planet (m/(s*s))

final PVector E = new PVector(0.0, 0.0, 0.0); //star

final PVector s0p = new PVector(150E9, 0.0, 0.0); // planet
final PVector v0p = new PVector(0, 30000, 0.0); // planet

// Main code:
// Converts distances from world length to pixel length
float worldToPixels(float dist)
{
  return dist*PIXELS_PER_METER;
}

// Converts distances from pixel length to world length
float pixelsToWorld(float dist)
{
  return dist/PIXELS_PER_METER;
}

// Converts a point from world coordinates to screen coordinates
void worldToScreen(PVector worldPos, PVector screenPos)
{
  screenPos.x = 0.5*DISPLAY_SIZE_X + (worldPos.x - DISPLAY_CENTER.x)*PIXELS_PER_METER;
  screenPos.y = 0.5*DISPLAY_SIZE_Y - (worldPos.y - DISPLAY_CENTER.y)*PIXELS_PER_METER;
}

// Converts a point from screen coordinates to world coordinates
void screenToWorld(PVector screenPos, PVector worldPos)
{
  worldPos.x = ((screenPos.x - 0.5*DISPLAY_SIZE_X)/PIXELS_PER_METER) + DISPLAY_CENTER.x;
  worldPos.y = ((0.5*DISPLAY_SIZE_Y - screenPos.y)/PIXELS_PER_METER) + DISPLAY_CENTER.y;
}

void drawStaticEnvironment()
{
  background(BACKGROUND_COLOR[0], BACKGROUND_COLOR[1], BACKGROUND_COLOR[2]);

  strokeWeight(1);

  // Star
  fill(255,255,0);
  PVector screenPos = new PVector();
  worldToScreen(E, screenPos);
  circle(screenPos.x, screenPos.y, worldToPixels(OBJECTS_SIZE));
  
  // Planet initial position
  fill(0,255,94);
  PVector screenPos2 = new PVector();
  worldToScreen(s0p, screenPos2);
  circle(screenPos2.x, screenPos2.y, worldToPixels(OBJECTS_SIZE));
  
}

void drawMovingElements()
{
  fill(0,145,255);
  strokeWeight(1);

  PVector screenPos = new PVector();
  worldToScreen(_s, screenPos);
  circle(screenPos.x, screenPos.y, worldToPixels(OBJECTS_SIZE));
}

void PrintInfo()
{
  println("Energy: " + _energy + " J");
  println("Elapsed time = " + _elapsedTime + " s");
  println("Simulated time = " + _simTime + " s \n");
}

void initSimulation()
{
  _simTime = 0.0;
  _elapsedTime = 0.0;
  
  _s = s0p.copy();
  _v = v0p.copy();
  _a.set(0.0,0.0,0.0);
}

void updateSimulation()
{
  switch (integrator)
  {
  case EXPLICIT_EULER:
    updateSimulationExplicitEuler();
    break;

  case SIMPLECTIC_EULER:
    updateSimulationSimplecticEuler();
    break;

  case HEUN:
    updateSimulationHeun();
    break;

  case RK2:
    updateSimulationRK2();
    break;

  case RK4:
    updateSimulationRK4();
    break;
  }
  
  _simTime += SIM_STEP;
  
  _output.println("Tiempo " + _simTime + "Energía " + _energy + "\n");

}

void updateSimulationExplicitEuler()
{
   _a = calculateAcceleration(_s, _v);
  
  _s.add(PVector.mult(_v, SIM_STEP));
  _v.add(PVector.mult(_a, SIM_STEP));
}

void updateSimulationSimplecticEuler()
{
  _a = calculateAcceleration(_s, _v);
  
  _v.add(PVector.mult(_a, SIM_STEP));
  _s.add(PVector.mult(_v, SIM_STEP));
}

void updateSimulationHeun()
{
  // Integración numérica de la velocidad
  // Calcular aceleración _a(s_i, v_i)
  _a = calculateAcceleration(_s, _v);
  // Paso de Euler -> actualizo s2, v2
  PVector s2 = PVector.add(_s, PVector.mult(_v, SIM_STEP));
  PVector v2 = PVector.add(_v, PVector.mult(_a, SIM_STEP));
  
  PVector v_promedio = PVector.mult(PVector.add(_v, v2),0.5);
  _s.add(PVector.mult(v_promedio, SIM_STEP));
  
  // Integrar la aceleracion
  // Calcular la aceleracon al final del intervalo
  // a2 = a(s2, v2)
  PVector a2 = calculateAcceleration(s2, v2);
  // Promedio de aceleraciones --> (_a, a2)
  PVector a_promedio = PVector.mult(PVector.add(_a, a2), 0.5);
  //Actualizar la velocidad _v, con la aceleracion promedia
  _v.add(PVector.mult(a_promedio, SIM_STEP));
}

void updateSimulationRK2()
{
  _a = calculateAcceleration(_s, _v);
  // Metodo original:
  // k1s = v(t)*h
  PVector k1s = PVector.mult(_v, SIM_STEP);
  // k1v = a(s(t),v(t))*h
  PVector k1v = PVector.mult(_a, SIM_STEP);
  
  PVector s2 = PVector.add(_s, PVector.mult(k1s, 0.5));
  PVector v2 = PVector.add(_v, PVector.mult(k1v, 0.5));
  PVector a2 = calculateAcceleration(s2,v2);
  // k2v = a(s(t)+k1s/2, v(t)+k1v/2)*h
  PVector k2v = PVector.mult(a2, SIM_STEP);
  // k2s = (v(t)+k1v/2)*h
  PVector k2s = PVector.mult(PVector.add(_v, PVector.mult(k1v, 0.5)), SIM_STEP);
  
  _v.add(k2v);
  _s.add(k2s);
}

void updateSimulationRK4()
{
  _a = calculateAcceleration(_s, _v);
  // k1v = a(s(t),v(t))*h
  PVector k1v = PVector.mult(_a, SIM_STEP);
  // k1s = v(t)*h
  PVector k1s = PVector.mult(_v, SIM_STEP);
  
  PVector s2 = PVector.add(_s, PVector.mult(k1s, 0.5));
  PVector v2 = PVector.add(_v, PVector.mult(k1v, 0.5));
  PVector a2 = calculateAcceleration(s2, v2);
  // k2v = a(s(t)+k1s/2, v(t)+k1v/2)*h
  PVector k2v = PVector.mult(a2, SIM_STEP);
  // k2s = (v(t)+k1v/2)*h
  PVector k2s = PVector.mult(PVector.add(_v, PVector.mult(k1v, 0.5)), SIM_STEP);

  PVector s3 = PVector.add(_s, PVector.mult(k2s, 0.5));
  PVector v3 = PVector.add (_v, PVector.mult(k2v, 0.5));
  PVector a3 = calculateAcceleration(s3,v3);
  // k3v = a(s(t)+k2s/2, v(t)+k2v/2)*h
  PVector k3v = PVector.mult(a3, SIM_STEP);
  // k3s = (v(t)+k2v/2)*h
  PVector k3s = PVector.mult(PVector.add(_v, PVector.mult(k2v, 0.5)), SIM_STEP);

  PVector s4 = PVector.add(_s, k3s);
  PVector v4 = PVector.add(_v, k3s);
  PVector a4 = calculateAcceleration(s4,v4);
  // k4v = a(s(t)+k3s, v(t)+k3v)*h
  PVector k4v = PVector.mult(a4, SIM_STEP);
  // k4s = (v(t)+k3v)*h
  PVector k4s = PVector.mult(PVector.add(_v,k3v), SIM_STEP);

  // v(t+h) = v(t) + (1/6)*k1v + (1/3)*k2v + (1/3)*k3v +(1/6)*k4v  
  // s(t+h) = s(t) + (1/6)*k1s + (1/3)*k2s + (1/3)*k3s +(1/6)*k4s  
  
  _v.add(PVector.mult(k1v, 1/6.0));
  _v.add(PVector.mult(k2v, 1/3.0));
  _v.add(PVector.mult(k3v, 1/3.0));
  _v.add(PVector.mult(k4v, 1/6.0));
  
  _s.add(PVector.mult(k1s, 1/6.0));
  _s.add(PVector.mult(k2s, 1/3.0));
  _s.add(PVector.mult(k3s, 1/3.0));
  _s.add(PVector.mult(k4s, 1/6.0));
}


PVector calculateAcceleration(PVector s, PVector v)
{ 
  PVector distancia = PVector.sub(E, s);
  float r = distancia.mag();
  float radio = r * r;
  float division = Gc*Me / radio;
  float formula = division * Mp;
  
  distancia.normalize();
  PVector Fg = PVector.mult(distancia, formula);

  PVector a = PVector.div(Fg, Mp);
  return a;
}


void calculateEnergy()
{  
  float Ek, Ep;
  float v = _v.mag();
  float g = G.mag();

  Ek = (Mp*v*v)*0.5;
  
  PVector distancia = PVector.sub(E, _s);
  float r = distancia.mag();
  
  float division = -g*Me / r;
  Ep = division * Mp;
  
  _energy = Ek + Ep;
}

void settings()
{
  if (FULL_SCREEN)
  {
    fullScreen();
    DISPLAY_SIZE_X = displayWidth;
    DISPLAY_SIZE_Y = displayHeight;
  } 
  else
    size(DISPLAY_SIZE_X, DISPLAY_SIZE_Y);
}

void setup()
{
  frameRate(DRAW_FREQ);
  _lastTimeDraw = millis();

  _output = createWriter(FILE_NAME);
  initSimulation();
}

void draw()
{
  int now = millis();
  _deltaTimeDraw = (now - _lastTimeDraw)/1000.0;
  _elapsedTime += _deltaTimeDraw;
  _lastTimeDraw = now;

  println("\nDraw step = " + _deltaTimeDraw + " s - " + 1.0/_deltaTimeDraw + " Hz");

  if (REAL_TIME)
  {
    float expectedSimulatedTime = TIME_AGES*_deltaTimeDraw;
    float expectedIterations = expectedSimulatedTime/SIM_STEP;
    int iterations = 0; 

    for (; iterations < floor(expectedIterations); iterations++)
      updateSimulation();

    if ((expectedIterations - iterations) > random(0.0, 1.0))
    {
      updateSimulation();
      iterations++;
    }

    println("Expected Simulated Time: " + expectedSimulatedTime);
    println("Expected Iterations: " + expectedIterations);
    println("Iterations: " + iterations);
  } 
  else
    updateSimulation();

  drawStaticEnvironment();
  drawMovingElements();

  text("Integrator: ", width-200, height-900);
  text(integrador, width-140, height-900);
  text("Paso de simulación: ", width-200, height-880);
  text(paso_sim, width-90, height-880);
  
  calculateEnergy();
  PrintInfo();
}

void mouseClicked() 
{
  PVector raton = new PVector(mouseX, mouseY, 0.0);
  PVector p = new PVector();
  
  screenToWorld(raton, p);  
  s0p.set(p.x, p.y, 0.0);
  
  initSimulation();
}

void keyPressed()
{
  if (key == 'e' || key == 'E')
  {
    integrator = IntegratorType.EXPLICIT_EULER; 
    updateSimulation();
    integrador = "Euler";
  }
  else if (key == 's' || key == 'S')
  {
    integrator = IntegratorType.SIMPLECTIC_EULER;
    updateSimulation();
    integrador = "Simplectic Euler";
  }
  else if (key == 'h' || key == 'H')
  {
    integrator = IntegratorType.HEUN;
    updateSimulation();
    integrador = "Heun";
  }
  else if (key == 'r' || key == 'R')
  {
    integrator = IntegratorType.RK2;
    updateSimulation();
    integrador = "RK2";
  }
  else if (key == 'k' || key == 'K')
  {
    integrator = IntegratorType.RK4;
    updateSimulation();
    integrador = "RK4";
  }
    
  else if (key == 'm' || key == 'M')
  {
    SIM_STEP = SIM_STEP + 0.1;
    paso_sim = str(SIM_STEP);
  }
    
  else if (key == 'n' || key == 'N')
  {
    if (SIM_STEP - 0.1 <= 0.1)
      SIM_STEP = 0.1;
    else
      SIM_STEP = SIM_STEP - 0.1;
    paso_sim = str(SIM_STEP);
  }
  // reset
  else if (key == 't' || key == 'T')
  {
    initSimulation();  
  }
    
  println("Key Pressend and change the integrator to: ", integrator);
}

void stop()
{
  _output.flush();
  _output.close();
}
