package coursework;

import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import model.Fitness;
import model.Individual;
import model.LunarParameters.DataSet;
import model.NeuralNetwork;

/**
 * Implements a basic Evolutionary Algorithm to train a Neural Network
 * 
 * You Can Use This Class to implement your EA or implement your own class that extends {@link NeuralNetwork} 
 * 
 */
public class ExampleEvolutionaryAlgorithm extends NeuralNetwork {
	

	/**
	 * The Main Evolutionary Loop
	 */
	@Override
	public void run() {		
		//Initialise a population of Individuals with random weights
		population = initialise();

		//Record a copy of the best Individual in the population
		best = getBest();
		System.out.println("Best From Initialisation " + best);

		/**
		 * main EA processing loop
		 */		
		
		while (evaluations < Parameters.maxEvaluations) {

			/**
			 * this is a skeleton EA - you need to add the methods.
			 * You can also change the EA if you want 
			 * You must set the best Individual at the end of a run
			 * 
			 */

			// Select 2 Individuals from the current population. Currently returns random Individual
			Individual parent1 = select(); 
			Individual parent2 = select();

			// Generate a child by crossover. Not Implemented			
			ArrayList<Individual> children = reproduce(parent1, parent2);			
			
			//mutate the offspring
			mutate(children);
			
			// Evaluate the children
			evaluateIndividuals(children);			

			// Replace children in population
			replace(children);

			// check to see if the best has improved
			best = getBest();
			
			// Implemented in NN class. 
			outputStats();
			
			//Increment number of completed generations			
		}

		//save the trained network to disk
		saveNeuralNetwork();
	}

	

	/**
	 * Sets the fitness of the individuals passed as parameters (whole population)
	 * 
	 */
	private void evaluateIndividuals(ArrayList<Individual> individuals) {
		for (Individual individual : individuals) {
			individual.fitness = Fitness.evaluate(individual, this);
		}
	}


	/**
	 * Returns a copy of the best individual in the population
	 * 
	 */
	private Individual getBest() {
		best = null;;
		for (Individual individual : population) {
			if (best == null) {
				best = individual.copy();
			} else if (individual.fitness < best.fitness) {
				best = individual.copy();
			}
		}
		return best;
	}

	/**
	 * Generates a randomly initialised population
	 * 
	 */
	private ArrayList<Individual> initialise() {
		population = new ArrayList<>();
		for (int i = 0; i < Parameters.popSize; ++i) {
			//chromosome weights are initialised randomly in the constructor
			Individual individual = new Individual();
			population.add(individual);
		}
		evaluateIndividuals(population);
		return population;
	}

	/**
	 * Selection --
	 * 
	 * Tournament Selection
	 */

	private Individual select() {	
		int i, picked, bestIndex;
		bestIndex=0;
		double bestFit=100;
		
		 for (i = 0; i < 10; i++)
         {
			picked = ThreadLocalRandom.current().nextInt(1,Parameters.popSize);
			if(population.get(picked).fitness <= bestFit)
			{
				bestFit = population.get(picked).fitness;
				bestIndex = picked;
			}
		}
			
		return population.get(bestIndex);
	}

	
 /*
  * 
  * Roulette Selection.
  * Referenced in report. Inferior results to Tournament selection
  * 
	private Individual select() {	
		double total_fitness = 0.0;
		ArrayList<Double> fitnessTable = new ArrayList<Double>();
		for(Individual i : population) {
			total_fitness += 1 - i.fitness;
			fitnessTable.add( 1 - i.fitness);
		}
		Random rd = new Random(); 
		double randomFitness = total_fitness  * rd.nextDouble();
		int idx = 0;
		double currentSum = 0;
		for(int i = 0; i < population.size(); i++) {
			currentSum += fitnessTable.get(i);
			if(currentSum > randomFitness) {
				idx = i;
				break;
			}
		}
		return population.get(idx); 
	}

*/


	/**
	 * Crossover / Reproduction
	 * 
	 * NEEDS REPLACED with proper method this code just returns exact copies of the
	 * parents. 
	 */
	

	private ArrayList<Individual> reproduce(Individual parent1, Individual parent2) {
		
		int crosspoint1,crosspoint2,j;
		ArrayList<Individual> children = new ArrayList<>();
		Individual child1= new Individual();
		Individual child2= new Individual();
		crosspoint1 = ThreadLocalRandom.current().nextInt(0,Parameters.getNumGenes() + 1);
		crosspoint2 = ThreadLocalRandom.current().nextInt(0,Parameters.getNumGenes() + 1);
		
		if (crosspoint1> crosspoint2){
			

			for (j = 0; j < crosspoint1; j++)
			{
				//Add the int at point j to the child
				child1.chromosome[j]=parent1.chromosome[j];
				child2.chromosome[j]=parent2.chromosome[j];
			}
			
			
			
			for (j = crosspoint1; j < crosspoint2; j++)
			{
				//Add the int at point j to the child
				child1.chromosome[j]=parent2.chromosome[j];
				child2.chromosome[j]=parent1.chromosome[j];
			}
			
			for (j = crosspoint2; j < Parameters.getNumGenes(); j++)
			{
				//Add the int at point j to the child
				child1.chromosome[j]=parent1.chromosome[j];
				child2.chromosome[j]=parent2.chromosome[j];
			}
	
			children.add(child1);
			children.add(child2);
			
		}
		
		else {
			

			for (j = 0; j < crosspoint2; j++)
			{
				//Add the int at point j to the child
				child1.chromosome[j]=parent1.chromosome[j];
				child2.chromosome[j]=parent2.chromosome[j];
			}
			
			
			
			for (j = crosspoint2; j < crosspoint1; j++)
			{
				//Add the int at point j to the child
				child1.chromosome[j]=parent2.chromosome[j];
				child2.chromosome[j]=parent1.chromosome[j];
			}
			
			for (j = crosspoint1; j < Parameters.getNumGenes(); j++)
			{
				//Add the int at point j to the child
				child1.chromosome[j]=parent1.chromosome[j];
				child2.chromosome[j]=parent2.chromosome[j];
			}
	
			children.add(child1);
			children.add(child2);
			
		}
		
		return children;
		
	} 

	
	/*One point crossover
		
	private ArrayList<Individual> reproduce(Individual parent1, Individual parent2) {
		
		int crosspoint,j;
		ArrayList<Individual> children = new ArrayList<>();
		Individual child1= new Individual();
		Individual child2= new Individual();
		crosspoint = ThreadLocalRandom.current().nextInt(0,Parameters.getNumGenes() + 1);
		for (j = 0; j < crosspoint; j++)
		{
			//Add the int at point j to the child
			child1.chromosome[j]=parent1.chromosome[j];
			child2.chromosome[j]=parent2.chromosome[j];
		}
		
		
		
		for (j = crosspoint; j < Parameters.getNumGenes(); j++)
		{
			//Add the int at point j to the child
			child1.chromosome[j]=parent2.chromosome[j];
			child2.chromosome[j]=parent1.chromosome[j];
		}

		children.add(child1);
		children.add(child2);

		return children;
	} 

	*/
	/**
	 * Mutation
	 * 
	 * 
	 */
	private void mutate(ArrayList<Individual> individuals) {		

		/*
		 * RSM Mutation.
		 * Referenced in report, inferior results to the standard mutation.
		 * 
		
		for(Individual individual : individuals) {
			Random random = new Random(); 
			int a = random.nextInt(individual.chromosome.length);
			Random random2 = new Random(); 
			int b = random2.nextInt(individual.chromosome.length - a) + a; 
			while(a < b) {
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[a] += (Parameters.mutateChange);
						individual.chromosome[b] += (Parameters.mutateChange);
					} else {
						individual.chromosome[a] -= (Parameters.mutateChange);
						individual.chromosome[b] -= (Parameters.mutateChange);
					}
					double temp = individual.chromosome[a];
					individual.chromosome[a] = individual.chromosome[b];
					individual.chromosome[b] = temp;
					Random random3 = new Random(); 
					int p = random3.nextInt(individual.chromosome.length); 
					if(p <= 0.5) {
						Random random4 = new Random();
						int j = random4.nextInt(individual.chromosome.length); 
						if (Parameters.random.nextBoolean()) {
							individual.chromosome[j] += (Parameters.mutateChange);

						} else {
							individual.chromosome[j] -= (Parameters.mutateChange);
						}
						double temp2 = individual.chromosome[a]; 
						individual.chromosome[a] = individual.chromosome[j];
						individual.chromosome[j] = temp2; 
					}


				}
				else {

				}
				a++;
				b--; 
			}
		}
	}
		
		*/
		/* Standard Mutation
		 * 
		 * 
		 *
		 * */
		for(Individual individual : individuals) {
			for (int i = 0; i < individual.chromosome.length; i++) {
				if (Parameters.random.nextDouble() < Parameters.mutateRate) {
					if (Parameters.random.nextBoolean()) {
						individual.chromosome[i] += (Parameters.mutateChange);
					} else {
						individual.chromosome[i] -= (Parameters.mutateChange);
					}
				}
			}
		}		
	}
	/*

	/**
	 * 
	 * Replaces the worst member of the population 
	 * (regardless of fitness)
	 * 
	 */
	private void replace(ArrayList<Individual> individuals) {
		Individual best = null;
		
		for(Individual individual : individuals) {
			if(best == null) {
				best = individual;
			}
			else if(individual.fitness < best.fitness) {
				best = individual;
			}
			
		}	
		int idx = getWorstIndex();		
		population.set(idx, best);
	}




	

	/**
	 * Returns the index of the worst member of the population
	 * @return
	 */
	private int getWorstIndex() {
		Individual worst = null;
		int idx = -1;
		for (int i = 0; i < population.size(); i++) {
			Individual individual = population.get(i);
			if (worst == null) {
				worst = individual;
				idx = i;
			} else if (individual.fitness > worst.fitness) {
				worst = individual;
				idx = i; 
			}
		}
		return idx;
	}	

	@Override
	public double activationFunction(double x) {
		if (x < -20.0) {
			return -1.0;
		} else if (x > 20.0) {
			return 1.0;
		}
		return Math.tanh(x);
	}
}
