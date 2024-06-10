package MLTKO;

import java.io.IOException;

public class testMLTKO {

	public static void main(String [] args) throws IOException {
		String	trans = "chainstore";
		String	dataset = trans + "_trans.txt";
		String	taxonomy = trans + "_tax.txt";
		int	K = 2000;
		//2436777
		//2441027
		boolean	eucp = true;
		boolean	bench =  false;					// benchmarking flag
		String	output = "output"+(eucp ? "_eucp" : "")+".txt";
		
		long	totalRun = 5;					// maximum number of evaluations per benchmark
		double	avgRuntime, minRuntime;			// for statistical purposes
		double	avgMemory, minMemory;
		
		// For testing with non-taxonomy dataset, it would reverted to TKO-Basic
		//String	taxonomy = "empty_tax.txt";

		if (!bench)								// regular run of the algorithm
		{
			//System.gc();
			AlgoMLTKO algo = new AlgoMLTKO(eucp);			
			algo.printLogo();
			algo.runAlgorithm(dataset, taxonomy, output, K);		
			algo.printStats();
		}
		else
		{
			avgRuntime = 0L; minRuntime = Integer.MAX_VALUE;
			avgMemory = 0L; minMemory = Double.MAX_VALUE;
			System.out.println("Benchmarking algorithm...");
			for (long i = 0; i < totalRun; i++)
			{
				System.gc();

				AlgoMLTKO algo = new AlgoMLTKO(eucp);
				algo.benchmark = true;
				
				System.out.print("Loop "+(i+1)+": ");
				Thread job = new Thread(new Runnable() {
				    @Override
				    public void run() {
				    	System.out.println("Job started.");
						try {
							algo.runAlgorithm(dataset, taxonomy, null, K);
						}
						catch (IOException e) { e.printStackTrace(); }
				    }
				});  
				job.setPriority(Thread.MAX_PRIORITY);
				job.start();
				try {
					job.join();
				} catch (InterruptedException e) { e.printStackTrace(); }
				System.out.println("Job finished.");
								
				algo.printStats();
				avgRuntime += algo.algoRuntime;
				avgMemory += algo.algoMemUsage;
				if (minRuntime > algo.algoRuntime) minRuntime = algo.algoRuntime;
				if (minMemory > algo.algoMemUsage) minMemory = algo.algoMemUsage;				
			}
			avgRuntime /= totalRun;
			avgMemory /= totalRun;
			
			System.out.println("===== BENCHMARK RESULTS =====");
			System.out.println("- Avg runtime : " + avgRuntime + " ms (" + avgRuntime / 1000.0 + " s)");
			System.out.println("- Avg memory  : " + avgMemory  + " MB");
			System.out.println("- Best runtime: " + minRuntime + " ms (" + minRuntime / 1000.0 + " s)");
			System.out.println("- Best memory : " + minMemory  + " MB");
			System.out.println("=============================");
		}
	}	
}
