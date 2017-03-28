//==============================================================================
/*
\author    Your Name
*/
//==============================================================================

//------------------------------------------------------------------------------
#include "chai3d.h"
//------------------------------------------------------------------------------
#include <GLFW/glfw3.h>
//------------------------------------------------------------------------------
using namespace chai3d;
using namespace std;
//------------------------------------------------------------------------------


struct Mass {
	cShapeSphere *p;
	double mass;		//the mass of particle (kg)
	cVector3d pos;		//position of particle	(x,y,z; meters)
	cVector3d vel;		//velocity of particle	(meters per second??)
	cVector3d f;		//force applied on particle	(Newtons)

	Mass(cShapeSphere *particle, double _mass, cVector3d position) {
		p = particle;
		mass = _mass;
		pos = position;
		vel = cVector3d(0, 0, 0);
		p->setLocalPos(pos);
	}

	void updateParticle() {
		p->setLocalPos(pos);
	}
};

struct Spring {
	cShapeLine *line;
	double k;
	Mass *m1;
	Mass *m2;

	double restLength;		//rest length
	double ksd;				//spring damping coefficient

	Spring(cShapeLine *spring_line, double stiffness_k, double rest_length, double damping_coefficient_ksd, Mass *m_1, Mass *m_2) {
		line = spring_line;
		k = stiffness_k;
		restLength = rest_length;
		ksd = damping_coefficient_ksd;
		m1 = m_1;
		m2 = m_2;

		line->m_pointA = m1->pos;
		line->m_pointB = m2->pos;
	}

	Mass* getMass1() {
		return m1;
	}

	Mass* getMass2() {
		return m2;
	}

	void updateLine() {
		line->m_pointA = m1->pos;
		line->m_pointB = m2->pos;
	}
};

std::vector<Mass*> pActive;
std::vector<Mass*> pStatic;
std::vector<Spring*> pSpring;

void updateForceParticles();
void setupScene(int i);

bool DEBUG = false;


//------------------------------------------------------------------------------
// GENERAL SETTINGS
//------------------------------------------------------------------------------

// stereo Mode
/*
C_STEREO_DISABLED:            Stereo is disabled
C_STEREO_ACTIVE:              Active stereo for OpenGL NVDIA QUADRO cards
C_STEREO_PASSIVE_LEFT_RIGHT:  Passive stereo where L/R images are rendered next to each other
C_STEREO_PASSIVE_TOP_BOTTOM:  Passive stereo where L/R images are rendered above each other
*/
cStereoMode stereoMode = C_STEREO_DISABLED;

// fullscreen mode
bool fullscreen = false;

// mirrored display
bool mirroredDisplay = false;


//------------------------------------------------------------------------------
// DECLARED VARIABLES
//------------------------------------------------------------------------------

// a world that contains all objects of the virtual environment
cWorld* world;

// a camera to render the world in the window display
cCamera* camera;

// a light source to illuminate the objects in the world
cDirectionalLight *light;

// a haptic device handler
cHapticDeviceHandler* handler;

// a pointer to the current haptic device
cGenericHapticDevicePtr hapticDevice;

// a label to display the rates [Hz] at which the simulation is running
cLabel* labelRates;

// a small sphere (cursor) representing the haptic device 
cShapeSphere* cursor;

// flag to indicate if the haptic simulation currently running
bool simulationRunning = false;

// flag to indicate if the haptic simulation has terminated
bool simulationFinished = false;

// a frequency counter to measure the simulation graphic rate
cFrequencyCounter freqCounterGraphics;

// a frequency counter to measure the simulation haptic rate
cFrequencyCounter freqCounterHaptics;

// haptic thread
cThread* hapticsThread;

// a handle to window display context
GLFWwindow* window = NULL;

// current width of window
int width = 0;

// current height of window
int height = 0;

// swap interval for the display context (vertical synchronization)
int swapInterval = 1;


//------------------------------------------------------------------------------
// DECLARED FUNCTIONS
//------------------------------------------------------------------------------

// callback when the window display is resized
void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height);

// callback when an error GLFW occurs
void errorCallback(int error, const char* a_description);

// callback when a key is pressed
void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods);

// this function renders the scene
void updateGraphics(void);

// this function contains the main haptics simulation loop
void updateHaptics(void);

// this function closes the application
void close(void);


//==============================================================================
/*
TEMPLATE:    application.cpp

Description of your application.
*/
//==============================================================================

int main(int argc, char* argv[])
{
	//--------------------------------------------------------------------------
	// INITIALIZATION
	//--------------------------------------------------------------------------

	cout << endl;
	cout << "-----------------------------------" << endl;
	cout << "CHAI3D" << endl;
	cout << "-----------------------------------" << endl << endl << endl;
	cout << "Keyboard Options:" << endl << endl;
	cout << "[f] - Enable/Disable full screen mode" << endl;
	cout << "[m] - Enable/Disable vertical mirroring" << endl;
	cout << "[q] - Exit application" << endl;
	cout << endl << endl;


	//--------------------------------------------------------------------------
	// OPENGL - WINDOW DISPLAY
	//--------------------------------------------------------------------------

	// initialize GLFW library
	if (!glfwInit())
	{
		cout << "failed initialization" << endl;
		cSleepMs(1000);
		return 1;
	}

	// set error callback
	glfwSetErrorCallback(errorCallback);

	// compute desired size of window
	const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
	int w = 0.8 * mode->height;
	int h = 0.5 * mode->height;
	int x = 0.5 * (mode->width - w);
	int y = 0.5 * (mode->height - h);

	// set OpenGL version
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);

	// set active stereo mode
	if (stereoMode == C_STEREO_ACTIVE)
	{
		glfwWindowHint(GLFW_STEREO, GL_TRUE);
	}
	else
	{
		glfwWindowHint(GLFW_STEREO, GL_FALSE);
	}

	// create display context
	window = glfwCreateWindow(w, h, "CHAI3D", NULL, NULL);
	if (!window)
	{
		cout << "failed to create window" << endl;
		cSleepMs(1000);
		glfwTerminate();
		return 1;
	}

	// get width and height of window
	glfwGetWindowSize(window, &width, &height);

	// set position of window
	glfwSetWindowPos(window, x, y);

	// set key callback
	glfwSetKeyCallback(window, keyCallback);

	// set resize callback
	glfwSetWindowSizeCallback(window, windowSizeCallback);

	// set current display context
	glfwMakeContextCurrent(window);

	// sets the swap interval for the current display context
	glfwSwapInterval(swapInterval);

#ifdef GLEW_VERSION
	// initialize GLEW library
	if (glewInit() != GLEW_OK)
	{
		cout << "failed to initialize GLEW library" << endl;
		glfwTerminate();
		return 1;
	}
#endif


	//--------------------------------------------------------------------------
	// WORLD - CAMERA - LIGHTING
	//--------------------------------------------------------------------------

	// create a new world.
	world = new cWorld();

	// set the background color of the environment
	world->m_backgroundColor.setBlack();

	// create a camera and insert it into the virtual world
	camera = new cCamera(world);
	world->addChild(camera);

	// position and orient the camera
	camera->set(cVector3d(0.2, 0.0, 0.0),    // camera position (eye)
		cVector3d(0.0, 0.0, 0.0),    // look at position (target)
		cVector3d(0.0, 0.0, 1.0));   // direction of the (up) vector

									 // set the near and far clipping planes of the camera
	camera->setClippingPlanes(0.01, 10.0);

	// set stereo mode
	camera->setStereoMode(stereoMode);

	// set stereo eye separation and focal length (applies only if stereo is enabled)
	camera->setStereoEyeSeparation(0.01);
	camera->setStereoFocalLength(0.5);

	// set vertical mirrored display mode
	camera->setMirrorVertical(mirroredDisplay);

	// create a directional light source
	light = new cDirectionalLight(world);

	// insert light source inside world
	world->addChild(light);

	// enable light source
	light->setEnabled(true);

	// define direction of light beam
	light->setDir(-1.0, 0.0, 0.0);

	// create a sphere (cursor) to represent the haptic device
	cursor = new cShapeSphere(0.01);

	// insert cursor inside world
	world->addChild(cursor);


	//--------------------------------------------------------------------------
	// HAPTIC DEVICE
	//--------------------------------------------------------------------------

	// create a haptic device handler
	handler = new cHapticDeviceHandler();

	// get a handle to the first haptic device
	handler->getDevice(hapticDevice, 0);

	// open a connection to haptic device
	hapticDevice->open();

	// calibrate device (if necessary)
	hapticDevice->calibrate();

	// retrieve information about the current haptic device
	cHapticDeviceInfo info = hapticDevice->getSpecifications();

	// display a reference frame if haptic device supports orientations
	if (info.m_sensedRotation == true)
	{
		// display reference frame
		cursor->setShowFrame(true);

		// set the size of the reference frame
		cursor->setFrameSize(0.05);
	}

	// if the device has a gripper, enable the gripper to simulate a user switch
	hapticDevice->setEnableGripperUserSwitch(true);

	// *** Setup scene 1
	setupScene(1);

	//--------------------------------------------------------------------------
	// WIDGETS
	//--------------------------------------------------------------------------

	// create a font
	cFontPtr font = NEW_CFONTCALIBRI20();

	// create a label to display the haptic and graphic rates of the simulation
	labelRates = new cLabel(font);
	labelRates->m_fontColor.setWhite();
	camera->m_frontLayer->addChild(labelRates);


	//--------------------------------------------------------------------------
	// START SIMULATION
	//--------------------------------------------------------------------------

	// create a thread which starts the main haptics rendering loop
	hapticsThread = new cThread();
	hapticsThread->start(updateHaptics, CTHREAD_PRIORITY_HAPTICS);

	// setup callback when application exits
	atexit(close);


	//--------------------------------------------------------------------------
	// MAIN GRAPHIC LOOP
	//--------------------------------------------------------------------------

	// call window size callback at initialization
	windowSizeCallback(window, width, height);

	// main graphic loop
	while (!glfwWindowShouldClose(window))
	{
		// get width and height of window
		glfwGetWindowSize(window, &width, &height);

		// render graphics
		updateGraphics();

		// swap buffers
		glfwSwapBuffers(window);

		// process events
		glfwPollEvents();

		// signal frequency counter
		freqCounterGraphics.signal(1);
	}

	// close window
	glfwDestroyWindow(window);

	// terminate GLFW library
	glfwTerminate();

	// exit
	return 0;
}

//------------------------------------------------------------------------------

void windowSizeCallback(GLFWwindow* a_window, int a_width, int a_height)
{
	// update window size
	width = a_width;
	height = a_height;
}

//------------------------------------------------------------------------------

void errorCallback(int a_error, const char* a_description)
{
	cout << "Error: " << a_description << endl;
}

//------------------------------------------------------------------------------

void keyCallback(GLFWwindow* a_window, int a_key, int a_scancode, int a_action, int a_mods)
{
	// filter calls that only include a key press
	if (a_action != GLFW_PRESS)
	{
		return;
	}

	// option - exit
	else if ((a_key == GLFW_KEY_ESCAPE) || (a_key == GLFW_KEY_Q))
	{
		glfwSetWindowShouldClose(a_window, GLFW_TRUE);
	}

	// option - toggle fullscreen
	else if (a_key == GLFW_KEY_F)
	{
		// toggle state variable
		fullscreen = !fullscreen;

		// get handle to monitor
		GLFWmonitor* monitor = glfwGetPrimaryMonitor();

		// get information about monitor
		const GLFWvidmode* mode = glfwGetVideoMode(monitor);

		// set fullscreen or window mode
		if (fullscreen)
		{
			glfwSetWindowMonitor(window, monitor, 0, 0, mode->width, mode->height, mode->refreshRate);
			glfwSwapInterval(swapInterval);
		}
		else
		{
			int w = 0.8 * mode->height;
			int h = 0.5 * mode->height;
			int x = 0.5 * (mode->width - w);
			int y = 0.5 * (mode->height - h);
			glfwSetWindowMonitor(window, NULL, x, y, w, h, mode->refreshRate);
			glfwSwapInterval(swapInterval);
		}
	}

	// option - toggle vertical mirroring
	else if (a_key == GLFW_KEY_M)
	{
		mirroredDisplay = !mirroredDisplay;
		camera->setMirrorVertical(mirroredDisplay);
	}
}

//------------------------------------------------------------------------------

void close(void)
{
	// stop the simulation
	simulationRunning = false;

	// wait for graphics and haptics loops to terminate
	while (!simulationFinished) { cSleepMs(100); }

	// close haptic device
	hapticDevice->close();

	// delete resources
	delete hapticsThread;
	delete world;
	delete handler;
}

//------------------------------------------------------------------------------

void updateGraphics(void)
{
	/////////////////////////////////////////////////////////////////////
	// UPDATE WIDGETS
	/////////////////////////////////////////////////////////////////////

	// update haptic and graphic rate data
	labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " +
		cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");

	// update position of label
	labelRates->setLocalPos((int)(0.5 * (width - labelRates->getWidth())), 15);


	/////////////////////////////////////////////////////////////////////
	// RENDER SCENE
	/////////////////////////////////////////////////////////////////////

	// update shadow maps (if any)
	world->updateShadowMaps(false, mirroredDisplay);

	// render world
	camera->renderView(width, height);

	// wait until all GL commands are completed
	glFinish();

	// check for any OpenGL errors
	GLenum err;
	err = glGetError();
	if (err != GL_NO_ERROR) cout << "Error:  %s\n" << gluErrorString(err);
}

//------------------------------------------------------------------------------

void updateHaptics(void)
{
	// simulation in now running
	simulationRunning = true;
	simulationFinished = false;

	//world timer
	cPrecisionClock timer;
	timer.start();

	// main haptic simulation loop
	while (simulationRunning)
	{
		/////////////////////////////////////////////////////////////////////
		// READ HAPTIC DEVICE
		/////////////////////////////////////////////////////////////////////

		// read position 
		cVector3d position;
		hapticDevice->getPosition(position);

		// read orientation 
		cMatrix3d rotation;
		hapticDevice->getRotation(rotation);

		// read user-switch status (button 0)
		bool button = false;
		hapticDevice->getUserSwitch(0, button);


		/////////////////////////////////////////////////////////////////////
		// UPDATE 3D CURSOR MODEL
		/////////////////////////////////////////////////////////////////////

		// update position and orienation of cursor
		cursor->setLocalPos(position);
		cursor->setLocalRot(rotation);

		/////////////////////////////////////////////////////////////////////
		// COMPUTE FORCES
		/////////////////////////////////////////////////////////////////////

		cVector3d force(0, 0, 0);
		cVector3d torque(0, 0, 0);
		double gripperForce = 0.0;

		double delta_t = timer.getCurrentTimeSeconds();
		timer.start(true);	//reset the clock

		updateForceParticles();

		for (int i = 0; i < pActive.size(); i++) {
			Mass* m = pActive[i];

			cVector3d F_p = m->f;
			cVector3d acc = F_p / m->mass;

			cVector3d vel = m->vel + delta_t * acc;
			cVector3d pos = m->pos + delta_t * vel;

			m->vel = vel;
			m->pos = pos;
			m->updateParticle();
		}
		for (int i = 0; i < pSpring.size(); i++) {
			Spring* s = pSpring[i];
			s->updateLine();
		}

		/////////////////////////////////////////////////////////////////////
		// APPLY FORCES
		/////////////////////////////////////////////////////////////////////

		// send computed force, torque, and gripper force to haptic device
		hapticDevice->setForceAndTorqueAndGripperForce(force, torque, gripperForce);

		// signal frequency counter
		freqCounterHaptics.signal(1);
	}

	// exit haptics thread
	simulationFinished = true;
}

//------------------------------------------------------------------------------

/*
z
|
|
|
|________y
/
/
/
x
*/

void DebugMessage(std::string type, std::string output) {
	if (DEBUG == true) {
		std::cout << type << ": " << output << std::endl;
	}

}

//====================Forces ======================================//

cVector3d calculateForceGravity(Mass *m) {
	cVector3d F_gravity = m->mass * cVector3d(0.0, 0.0, -9.81);
	return F_gravity;
}

cVector3d calculateForceCollision(Mass *m) {
	// *** Need to fill it in
	// *** I don't think the particles should ever collide with other particles
	//so this will just be cursor-particle collision
	cVector3d F_collision(0.0, 0.0, 0.0);
	return F_collision;
}

cVector3d calculateForceDamping(Mass *m) {
	cVector3d F_damping;
	double c_air = 1.0;		//N/m
	F_damping = -c_air * m->vel;
	return F_damping;
}

cVector3d getForceOnParticle(Mass *m) {
	cVector3d Fg = calculateForceGravity(m);	//gravity
	cVector3d Fc = calculateForceCollision(m);	//collision
	cVector3d Fd = calculateForceDamping(m);	//air damping

												//cVector3d F = Fg + Fc + Fd;		//*** Add gravity, or take it out?
	cVector3d F = Fc + Fd;
	return F;
}

cVector3d calculateForceSpring(Spring *s) {
	Mass* m1 = s->getMass1();
	Mass* m2 = s->getMass2();

	double curr_length = (m1->pos - m2->pos).length();
	double rest_length = s->restLength;
	cVector3d dir = cNormalize(m1->pos - m2->pos);
	cVector3d F_spring_12 = -s->k * (curr_length - rest_length) * dir;

	return F_spring_12;
}

cVector3d calculateForceSpringDamping(Spring *s) {
	//*** need to check for correctness for springs
	Mass* m1 = s->getMass1();
	Mass* m2 = s->getMass2();


	double ksd = s->ksd;
	cVector3d vel_12 = m1->vel - m2->vel;
	cVector3d F_dir = cNormalize(m1->f);
	if (F_dir * F_dir > 0) {
		cVector3d F_springDamp = -ksd * ((vel_12 * F_dir) / (F_dir * F_dir)) * F_dir;
		return F_springDamp;
	}
	else {
		return cVector3d(0.0, 0.0, 0.0);
	}
}

cVector3d getForceFromSpring(Spring *s) {
	cVector3d Fs = calculateForceSpring(s);
	cVector3d Fsd = calculateForceSpringDamping(s);

	cVector3d F = Fs + Fsd;
	DebugMessage("FS", Fs.str());
	DebugMessage("Fsd", Fsd.str());

	return F;
}

void updateForceParticles() {
	//Update forces on particle
	for (int i = 0; i < pActive.size(); i++) {
		Mass* m = pActive[i];
		cVector3d F_particle = getForceOnParticle(m);
		//DebugMessage("F_PARTICLE", F_particle.str());
		m->f = F_particle;
	}


	for (int i = 0; i < pSpring.size(); i++) {
		Spring* s = pSpring[i];
		cVector3d F_spring = getForceFromSpring(s);
		DebugMessage("F_SPRING", F_spring.str());

		Mass* m1 = s->getMass1();
		Mass* m2 = s->getMass2();
		m1->f = m1->f + F_spring;
		m2->f = m2->f - F_spring;
	}
}

//======================= Scenes ======================================//

void addParticles(int size, double length, double radius, double mass) {

	cVector3d start_pos = cVector3d(0.0, -length / 2, 0.0);

	//active particles
	for (int i = 0; i < size; i++) {
		cShapeSphere* p = new cShapeSphere(radius);
		cVector3d interval = cVector3d(0.0, ((double)i + 1) * length / (size + 1), 0.0);
		Mass* m = new Mass(p, mass, start_pos + interval);
		p->m_material->setBlueLightSteel();
		world->addChild(p);
		pActive.push_back(m);
	}

	//static particles
	for (int i = 0; i < 2; i++) {
		cShapeSphere* p = new cShapeSphere(radius);
		cVector3d interval = cVector3d(0.0, (double)i *length, 0.0);
		Mass* m = new Mass(p, mass, start_pos + interval);
		p->m_material->setRedDark();
		world->addChild(p);
		pStatic.push_back(m);
	}
}

void addSprings(double k, double rest_length, double ksd) {

	for (int i = 0; i < pActive.size() - 1; i++) {
		Mass* m1 = pActive[i];
		Mass* m2 = pActive[i + 1];
		cShapeLine* l = new cShapeLine(m1->pos, m2->pos);

		Spring* s = new Spring(l, k, rest_length, ksd, m1, m2);

		l->m_material->setGreenLight();
		world->addChild(l);
		pSpring.push_back(s);
	}


	//add the two springs at the end
	for (int i = 0; i < 2; i++) {
		Mass* m1 = pStatic[i];
		Mass* m2 = pActive[i * (pActive.size() - 1)];
		cShapeLine* l = new cShapeLine();

		Spring* s = new Spring(l, k, rest_length, ksd, m1, m2);
		l->m_material->setGreenLight();
		world->addChild(l);
		pSpring.push_back(s);
	}
}

void setupScene1() {
	//mass
	int size = 1;
	double length = 0.1;
	double radius = 0.01;
	double mass = 1.0;
	addParticles(size, length, radius, mass);


	// springs
	double k = 400.0;
	double rest_length = 0.01;
	double ksd = 10.0;

	addSprings(k, rest_length, ksd);
}

void setupScene(int i) {
	if (i == 1) { setupScene1(); }
}

//===========================Cursor ====================================//
