//==============================================================================
/*
\author    Philmo Gu, Jackie Luc
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
	double k;			//contact stiffness

	Mass(cShapeSphere *particle, double _mass, cVector3d position) {
		p = particle;
		mass = _mass;
		pos = position;
		vel = cVector3d(0, 0, 0);
		k = 400;
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

std::vector<double> E_string = {330, 349, 370, 392, 415, 466, 494, 523, 554, 587, 622, 659, 698, 740, 784, 831, 880, 932, 988, 1047 };
std::vector<double> B_string = {247, 262, 277, 294, 311, 330, 349, 370, 392, 415};
std::vector<double> G_string = {196, 208, 220, 233, 247, 262, 277, 294, 311, 330, 349, 370};
std::vector<double> D_string = {147, 156, 165, 175, 185, 196, 208, 220, 233, 247, 262};



void updateForceParticles(cVector3d cursorPos);
void setupScene(int i);
void setupScene1();

bool DEBUG = false;

// background for the scene
cBackground* background;

// enable/disable widget to display graphics/haptics rate
bool showWidget = true;

//----sounds ----
//audio device to play sound
cAudioDevice* audioDevice;

//audio buffer to store sound files
std::vector<cAudioBuffer*> audioBuffer;
std::vector<double>	audioBufferFrequency;		//for tuning beyond the integer place

//position where the audio comes from (1 per string)
std::vector<cAudioSource*> audioSource;

void initSound();
void updateAudioSourcePosVel(cAudioSource* audioSource, Mass* m);
void changeFrequency(double delta_f, int buffer_index);
void changeNote(int note_pos, int buffer_index);

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

// a multi-mesh representing the haptic device as a guitar pick
cMultiMesh* cursor;

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
	cout << "[p] - Enable/Disable full screen mode" << endl;
	cout << "[m] - Enable/Disable vertical mirroring" << endl;
	cout << "[0] - Enable/Disable debug of display/haptic rate" << endl;
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

	// a multi-mesh representing the haptic device as a guitar pick
	cursor = new cMultiMesh();
	cursor->loadFromFile("GuitarPick.3ds");
	cursor->scale(0.01);

	// insert cursor inside world
	world->addChild(cursor);

	//==========================================================================
	// Inspiration from 19-space.cpp under chai3d-3.2.0/examples/GLFW/19-space
	//==========================================================================

	// create a background
	background = new cBackground();
	camera->m_backLayer->addChild(background);

	// load a texture file
	bool fileload = background->loadFromFile("images/guitar.jpg");
	if (!fileload)
	{
		#if defined(_MSVC)
		fileload = background->loadFromFile("images/guitar.jpg");
		#endif
	}
	if (!fileload)
	{
		cout << "Error - Background image failed to load correctly." << endl;
		close();
		return (-1);
	}

	//==========================================================================
	// End 19-space
	//==========================================================================
	
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
	//	SETUP AUDIO
	//	Code from: "25-sounds.cpp" of chai3d examples
	//--------------------------------------------------------------------------

	//create an audio buffer and load audio wave file
	initSound();

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

	double delta_f = 0.001;
	// filter calls that only include a key press
	if (a_action != GLFW_PRESS)
	{
		return;
	}

	// option - exit
	else if ((a_key == GLFW_KEY_ESCAPE))
	{
		glfwSetWindowShouldClose(a_window, GLFW_TRUE);
	}

	// option - toggle fullscreen
	else if (a_key == GLFW_KEY_P)
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
	//else if (a_key == GLFW_KEY_M)
	//{
	//	mirroredDisplay = !mirroredDisplay;
	//	camera->setMirrorVertical(mirroredDisplay);
	//}

	// reset scene
	//else if (a_key == GLFW_KEY_R)
	//{
	//	setupScene1();
	//}

	// toggle widget to enable/disable display of graphic / haptic rates
	else if (a_key == GLFW_KEY_0) { showWidget = !showWidget;}

	//tuning of guitar
	else if (a_key == GLFW_KEY_M) { changeFrequency(delta_f, 0); }
	else if (a_key == GLFW_KEY_N) { changeFrequency(-delta_f, 0); }

	else if (a_key == GLFW_KEY_J) { changeFrequency(delta_f, 1); }
	else if (a_key == GLFW_KEY_K) { changeFrequency(-delta_f, 1); }

	else if (a_key == GLFW_KEY_U) { changeFrequency(delta_f, 2); }
	else if (a_key == GLFW_KEY_I) { changeFrequency(-delta_f, 2); }

	else if (a_key == GLFW_KEY_7) { changeFrequency(delta_f, 3); }
	else if (a_key == GLFW_KEY_8) { changeFrequency(-delta_f, 3); }

	//changing notes of a guitar
	else if (a_key == GLFW_KEY_B) { changeNote(4, 0); }
	else if (a_key == GLFW_KEY_V) { changeNote(3, 0); }
	else if (a_key == GLFW_KEY_C) { changeNote(2, 0); }
	else if (a_key == GLFW_KEY_X) { changeNote(1, 0); }
	else if (a_key == GLFW_KEY_Z) { changeNote(0, 0); }
	
	else if (a_key == GLFW_KEY_G) { changeNote(4, 1); }
	else if (a_key == GLFW_KEY_F) { changeNote(3, 1); }
	else if (a_key == GLFW_KEY_D) { changeNote(2, 1); }
	else if (a_key == GLFW_KEY_S) { changeNote(1, 1); }
	else if (a_key == GLFW_KEY_A) { changeNote(0, 1); }

	else if (a_key == GLFW_KEY_T) { changeNote(4, 2); }
	else if (a_key == GLFW_KEY_R) { changeNote(3, 2); }
	else if (a_key == GLFW_KEY_E) { changeNote(2, 2); }
	else if (a_key == GLFW_KEY_W) { changeNote(1, 2); }
	else if (a_key == GLFW_KEY_Q) { changeNote(0, 2); }

	else if (a_key == GLFW_KEY_5) { changeNote(4, 3); }
	else if (a_key == GLFW_KEY_4) { changeNote(3, 3); }
	else if (a_key == GLFW_KEY_3) { changeNote(2, 3); }
	else if (a_key == GLFW_KEY_2) { changeNote(1, 3); }
	else if (a_key == GLFW_KEY_1) { changeNote(0, 3); }
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
	if (showWidget)
	{
		labelRates->setText(cStr(freqCounterGraphics.getFrequency(), 0) + " Hz / " +
			cStr(freqCounterHaptics.getFrequency(), 0) + " Hz");
	}
	else
		labelRates->setText("");


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

	int counter = 50000;
	int counterMin = 0;
	int counterMax = 100000;
	int counterDirection = 1;
	float prevVol = 0.0;

	// main haptic simulation loop
	while (simulationRunning)
	{
		/////////////////////////////////////////////////////////////////////
		// READ HAPTIC DEVICE
		/////////////////////////////////////////////////////////////////////

		// read position 
		cVector3d cursorPos;
		hapticDevice->getPosition(cursorPos);

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
		cursor->setLocalPos(cursorPos);
		cursor->setLocalRot(rotation);

		/////////////////////////////////////////////////////////////////////
		// COMPUTE FORCES
		/////////////////////////////////////////////////////////////////////

		//reset values
		cVector3d force(0, 0, 0);
		cVector3d torque(0, 0, 0);
		double gripperForce = 0.0;
		
		force = cVector3d(0, 0, 0);
			
		double delta_t = timer.getCurrentTimeSeconds();
		timer.start(true);	//reset the clock
					
		double* collisionDist = new double[4]();

		//if cursor position is near the mass spring
		for (int i = 0; i < pActive.size(); i++) {
			cVector3d m_pos = cVector3d(pActive[i]->pos.x(), 0.0, pActive[i]->pos.z());
			cVector3d c_pos = cVector3d(cursorPos.x(), 0.0, cursorPos.z());

			// get the distance between the mass position and the cursor position
			collisionDist[i] = (m_pos - c_pos).length();

			// we set the position of the mass in the same y position as the cursor
			// this allows us to compute collision forces quite easily
			if (collisionDist[i] < 0.01) {
				pActive[i]->pos = cVector3d(0.0, cursorPos.y(), pActive[i]->pos.z());
			}
		}

		updateForceParticles(cursorPos);
		
		// temporary virtual wall to find the strings easier
		if (cursorPos.x() < 0.01) {
			double Fx = -2000 * (cursorPos.x() - 0.01);
			force.add(cVector3d(Fx, 0.0, 0.0));
		}

		for (int i = 0; i < pActive.size(); i++) {
			Mass* m = pActive[i];

			cVector3d F_p = m->f;
			cVector3d acc = F_p / m->mass;

			cVector3d vel = m->vel + delta_t * acc;
			cVector3d pos = m->pos + delta_t * vel;

			// if there is collision, render the force to haptic device
			if (collisionDist[i] < 0.01) {
				force.add(-m->f);
				vel = cVector3d(0, 0, 0);
			}

			m->vel = vel;
			m->pos = pos;
			m->updateParticle();

			//----- audio (requires updated mass's velocity and position)

			if (collisionDist[i] > 0.01) {
					
				updateAudioSourcePosVel(audioSource[i], m);
					
				float maxAccleration = 100;
				float vol = std::min(abs(acc.z() / maxAccleration), 1.0);
				vol = std::max(vol, (float)0.001);
				if (vol != prevVol) {
					prevVol = vol;
					audioSource[i]->setGain(vol);
				}					
			}

			else {
				//minimum volume
				audioSource[i]->setGain(0.001);
			}

			
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

cVector3d calculateForceCollision(Mass *m, cVector3d cursorPos) {
	cVector3d F_collision(0, 0, 0);

	double dist = (m->pos - cursorPos).length();

	double collisionRadius = 0.01;
	double const scalar = 1.5;

	if (dist < collisionRadius) {
		double mag = scalar * m->k * (collisionRadius - dist);
		cVector3d d_vector = m->pos - cursorPos;

		// we are only concerned with the particle's behaviour in the y and z axis
		cVector3d dir = cNormalize(cVector3d(0.0, d_vector.y(), d_vector.z()));
		F_collision = mag * dir;
	}

	return F_collision;
}

cVector3d calculateForceDamping(Mass *m) {
	cVector3d F_damping;
	double c_air = 0.05;		//N/m
	F_damping = -c_air * m->vel;
	return F_damping;
}

cVector3d getForceOnParticle(Mass *m, cVector3d cursorPos) {
	cVector3d Fc = calculateForceCollision(m, cursorPos);
	cVector3d Fd = calculateForceDamping(m);	//air damping

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

	F_spring_12 = cVector3d(0, 0, F_spring_12.z());

	return F_spring_12;
}

cVector3d calculateForceSpringDamping(Spring *s) {
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
	return F;
}

void updateForceParticles(cVector3d cursorPos) {
	//Update forces on particle
	for (int i = 0; i < pActive.size(); i++) {
		Mass* m = pActive[i];
		cVector3d F_particle = getForceOnParticle(m, cursorPos);
		m->f = F_particle;
	}

	for (int i = 0; i < pSpring.size(); i++) {
		Spring* s = pSpring[i];
		cVector3d F_spring = getForceFromSpring(s);

		Mass* m1 = s->getMass1();
		Mass* m2 = s->getMass2();
		m1->f = m1->f + F_spring;
		m2->f = m2->f - F_spring;
	}
}

//======================= Scenes ======================================//

void addParticles(int size, double length, double radius, double mass) {


	cVector3d start_pos = cVector3d(0.0, (-length / 2) - 0.5, -0.015);

	for (int j = 0; j < size; j++) {

		//active particles
		cShapeSphere* p = new cShapeSphere(radius);
		cVector3d interval = cVector3d(0.0, (length / 2) + 0.5, 0.0);

		Mass* m = new Mass(p, mass, start_pos + interval);

		pActive.push_back(m);

		//static particles
		for (int i = 0; i < 2; i++) {
			cShapeSphere* p = new cShapeSphere(0.0015);
			cVector3d interval = cVector3d(0.0, (double)i *(length + 0.5225), 0.0);

			Mass* m = new Mass(p, mass, start_pos + interval);
			p->m_material->setWhite();

			// only graphically display the particle connected to the right-end of the guitar
			if (i == 1)
				world->addChild(p);

			pStatic.push_back(m);
		}

		start_pos += cVector3d(0.0, 0.0, 0.015);
	}
}

void addSprings(int size, double k, double rest_length, double ksd) {

	//add the two springs at the end
	for (int j = 0; j < size; j++) {
		for (int i = 0; i < 2; i++) {
			Mass* m1 = pStatic[2 * j + i];
			Mass* m2 = pActive[j];

			cShapeLine* l = new cShapeLine();

			Spring* s = new Spring(l, k, rest_length, ksd, m1, m2);
			l->m_material->setGreenLight();
			world->addChild(l);
			pSpring.push_back(s);
		}
	}
}

void clearScene() {
	world->clearAllChildren();
	world->addChild(camera);
	world->addChild(cursor);
	world->addChild(light);

	pActive.clear();
	pSpring.clear();
	pStatic.clear();
}

void setupScene1() {
	// clear the scene
	clearScene();

	//mass
	int size = 4;
	double length = 0.1;
	double radius = 0.01;
	double mass = 0.03;
	addParticles(size, length, radius, mass);


	// springs
	double k = 200.0;
	double rest_length = 0.01;
	double ksd = 0.0;
	addSprings(size, k, rest_length, ksd);
}

void setupScene(int i) {
	if (i == 1) { setupScene1(); }
}

//===========================Cursor ====================================//
/**	Code is now obselete. 
	Was intended to simulate plucking with the button on the Falcon.
	
int findNearestP(cVector3d cursorPos) {
	//find the nearest particle near the cursor
	int smallestIndex = 0;
	double smallestLength = 1000000.0;

	for (int i = 0; i < pActive.size(); i++) {
		double currentLength = (pActive[i]->pos - cursorPos).length();
		if (currentLength < smallestLength) {
			smallestLength = currentLength;
			smallestIndex = i;
		}
	}

	return smallestIndex;	//return the nearest index of pActive
}
*/

//========================Position =====================================//
/** Code is now obselete.
	Was intended to update the rest length after interactions.	

void updateRestLength(Spring* s1, Spring* s2, cVector3d cursorPos) {
	//get relative position of rest length based on cursor position
	double rest_length1 = fabs(s1->getMass1()->pos.y() - cursorPos.y()) / 5;
	double rest_length2 = fabs(s2->getMass1()->pos.y() - cursorPos.y()) / 5;

	s1->restLength = rest_length1;
	s2->restLength = rest_length2;
}
*/
//===========================Sound ====================================//

void initLoadedAudio(cAudioBuffer* audioBuffer, string fileName) {
	bool fileload = audioBuffer->loadFromFile(fileName);

	//check for errors
	if (!(fileload)) {
		cout << "Error - Sound file failed to load or initialize correctly." << endl;
		close();
	}
}

void createSinWave(double buf_size, double freq, double sample_rate, cAudioBuffer* audioBuffer) {
	double ampli = freq * 1500.0;		//high = good sounding high pitch tones, bad sounding low pitch tones
										//low = bad sounding high pitch tones, good sounding low pitch tones
	unsigned char *samples;
	samples = new unsigned char[buf_size];
	for (int i = 0; i < buf_size; ++i) {
		double a = ampli * sin((2.f*float(M_PI)*freq) / sample_rate * (double)i);
		samples[i] = a;
	}
	//We set up how fast the sound will be played (approximately) on here, eh?
	audioBuffer->setup(samples, buf_size, freq, false, sample_rate);
}

void createSquareWave(double buf_size, double freq, double sample_rate, cAudioBuffer* audioBuffer) {

	double ampli = 1.991;		//***somehow, second string goes missing???
	unsigned char *samples;
	samples = new unsigned char[buf_size];
	for (int i = 0; i < buf_size; ++i) {
		double a = sin((2.f*float(M_PI)*freq) / sample_rate * (double)i);
		if (a > 0) { samples[i] = ampli; }
		else { samples[i] = -ampli; }
		//samples[i] = a;
	}
	//We set up how fast the sound will be played (approximately) on here, eh?
	audioBuffer->setup(samples, buf_size, freq, false, sample_rate);
}

void initGeneratedAudio(cAudioBuffer* audioBuffer, float freq, unsigned sample_rate, int type) {
	type = 0;

	/* Fill buffer with Sine-Wave */
	int seconds = 4;
	size_t buf_size = seconds * sample_rate;
	if (type == 0) {
		createSinWave(buf_size, freq, sample_rate, audioBuffer);
	}
	if (type == 1) {
		createSquareWave(buf_size, freq, sample_rate, audioBuffer);
	}
}


void playAudioSource(cAudioSource* audioSource, cAudioBuffer* audioBuffer, cVector3d position, bool loop) {
	//assign audio buffer to audio source
	audioSource->setAudioBuffer(audioBuffer);

	//loop sound
	audioSource->setLoop(loop);

	//volume of sound
	audioSource->setGain(0.0);

	//pitch of sound
	audioSource->setPitch(1.0);

	//play sound
	audioSource->play();
}

void updateAudioSourcePosVel(cAudioSource* audioSource, Mass* m) {
	audioSource->setSourcePos(m->pos);
	audioSource->setSourceVel(m->vel);
}


void initSound() {
	//audio device to play sounds
	audioDevice = new cAudioDevice();

	//attach audio device to camera?
	camera->attachAudioDevice(audioDevice);

	//If it sounds like a guitar, tuned like a guitar, feels like a guitar,
	//but doesn't look like a guitar, is it a guitar?
	audioBufferFrequency.push_back(329.63);
	audioBufferFrequency.push_back(246.94);
	audioBufferFrequency.push_back(196.06);
	audioBufferFrequency.push_back(146.83);
		
	int size = 4;
	for (int i = 0; i < size; i++) {
		audioBuffer.push_back(new cAudioBuffer());
		initGeneratedAudio(audioBuffer[i], audioBufferFrequency[i], 1, 0);	//good (kinda static-y)
	}
	
	//create audio source (ie. where the audio is played from)
	for (int i = 0; i < size; i++) {
		audioSource.push_back(new cAudioSource());
		playAudioSource(audioSource[i], audioBuffer[i], cursor->getLocalPos(), true);
	}

}

void changeFrequency(double delta_f, int buffer_index) {

	audioBufferFrequency[buffer_index] = audioBufferFrequency[buffer_index] + delta_f;
	cout << "Freq[" << buffer_index << "]:" << audioBufferFrequency[buffer_index] << endl;
	audioSource[buffer_index]->~cAudioSource();
	audioSource[buffer_index] = new cAudioSource();
	initGeneratedAudio(audioBuffer[buffer_index], audioBufferFrequency[buffer_index], 1, 0);
	playAudioSource(audioSource[buffer_index], audioBuffer[buffer_index], cursor->getLocalPos(), true);
	audioSource[buffer_index]->play();
}

void changeNote(int note_pos, int buffer_index) {
	double freq;
	if (buffer_index == 0) { freq = E_string[note_pos] / audioBufferFrequency[buffer_index]; }
	else if (buffer_index == 1) { freq = B_string[note_pos] / audioBufferFrequency[buffer_index]; }
	else if (buffer_index == 2) { freq = G_string[note_pos] / audioBufferFrequency[buffer_index]; }
	else if (buffer_index == 3) { freq = D_string[note_pos] / audioBufferFrequency[buffer_index]; }
	cout << "Freq[" << note_pos << "," << buffer_index << "]: " << freq << endl;
	audioSource[buffer_index]->setPitch(freq);
}
