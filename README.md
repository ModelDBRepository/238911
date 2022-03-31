README for the model used in:

Barth, B. B., Henriquez, C. S., Grill, W. M., and Shen, X. (2017). "Electrical stimulation of gut motility guided by an in silico model" J Neural Eng 14(6).

ABSTRACT:
OBJECTIVE: Neuromodulation of the central and peripheral nervous systems is becoming increasingly important for treating a diverse set of diseases-ranging from Parkinson's Disease and epilepsy to chronic pain. However, neuromodulation of the gastrointestinal (GI) tract has achieved relatively limited success in treating functional GI disorders, which affect a significant population, because the effects of stimulation on the enteric nervous system (ENS) and gut motility are not well understood. Here we develop an integrated neuromechanical model of the ENS and assess neurostimulation strategies for enhancing gut motility, validated by in vivo experiments. APPROACH: The computational model included a network of enteric neurons, smooth muscle fibers, and interstitial cells of Cajal, which regulated propulsion of a virtual pellet in a model of gut motility. MAIN RESULTS: Simulated extracellular stimulation of ENS-mediated motility revealed that sinusoidal current at 0.5 Hz was more effective at increasing intrinsic peristalsis and reducing colon transit time than conventional higher frequency rectangular current pulses, as commonly used for neuromodulation therapy. Further analysis of the model revealed that the 0.5 Hz sinusoidal currents were more effective at modulating the pacemaker frequency of interstitial cells of Cajal. To test the predictions of the model, we conducted in vivo electrical stimulation of the distal colon while measuring bead propulsion in awake rats. Experimental results confirmed that 0.5 Hz sinusoidal currents were more effective than higher frequency pulses at enhancing gut motility. SIGNIFICANCE: This work demonstrates an in silico GI neuromuscular model to enable GI neuromodulation parameter optimization and suggests that low frequency sinusoidal currents may improve the efficacy of GI pacing.

This work was submitted by Bradley Barth (Duke University), and it relies on the Brian 2 simulator available in Python 2.7.
Brian 2 [Brian 2.0rc3] is available at https://github.com/brian-team/brian2/releases/tag/2.0rc3. 


INSTRUCTIONS:
First, install Python 2.7 and Brian 2. Instructions for installing Python 2.7 and Brian 2 can be found at https://www.python.org/ and https://github.com/brian-team/brian2/releases/tag/2.0rc3, respectively. Second, download AH_Neuron.py and Motility_Model.py. Third, run either python file in preferred IDE such as Spyder.

The AH_Neuron.py file demonstrates a neuron with and without the afterhyperpolarization current used in the model. It produces two figures: Cell-1_without-AH-Current.png and Cell-2_with-AH-Current.png, which are included under /Figures/ and reproduce data represented in Supplementary Figure 9a of Barth et al. 2019.

The Motility_Model.py file contains five modules: ExtracellStim(), pellet(), deterministic(), stochastic(), and Motility Model_CEs(). Running the model will output: (1) a summary figure, similar to Fig. 1c, (2) the transit time, and (3) the waveform used to stimulate the model. If the transit time is "Incomplete", the pellet did not reach the end of the model within the total run time. 


NOTES:           
To run the base model, run the Motility_Model.py file, then call the following in the command line:
	controlStim = ExtracellStim("Control")
	control = MotilityModel_CES(controlStim,"moderate")

ExtracellStim() is a module to create a stimulation vector based on the stimulation pattern and frequency (among other parameters).

Calling "pellet()" in the command line will reproduce the data represented in Supplementary Figure 9c of Barth et al. 2019, included under /Figures/Pellet-velocity_vs_pellet-size.png.

The model comparison between sine wave and pulsatile stimulation can be simulated by calling either "deterministic()" or "stochastic()". The difference between these two modules is the absence or presence of spontaneous junction potentials. Calling "deterministic()" will reproduce the data represented in Figure 3d and 3e of Barth et al. 2019, included under /Figures/deterministic-model.png. Calling "stochastic("moderate")" (identical to calling "stochastic()") will reproduce the data represented in Figure 5a and 5b of Barth et al. 2019, included under /Figures/stochastic-moderate-model.png. Calling "stochastic("high")" will reproduce the data represented in Figure 5c and 5d of Barth et al. 2019, included under /Figures/stochastic-high-model.png.
 
 
 
ADVANCED SETTINGS: 
The input parameter "advanced" can be "low", "moderate", or "high", and it refers to additional settings in the model described in Figure 5.

The base model uses has the "advanced" parameter setting "low". The "moderate" setting includes additional neural and muscular mechanims: increased junction potential amplitudes, ascending inhibition and descending excitation pathways, fiber conduction delays, after-hyperpolarization current and graded stretch response in intrinsic sensory neurons, and stochastic Poisson events to model spontaneous junction potentials.

The "high" advanced setting includes the same mechanisms as in the "moderate" advanced settings, except the synaptic weights of the additional pathways are increased by 60%.
