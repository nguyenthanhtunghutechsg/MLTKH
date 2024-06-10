package MLTKH;

import java.io.IOException;

public class testMLTKH {

	public static void main(String [] args) throws IOException {
		String	trans = "chainstore";
		String	dataset = trans + "_trans.txt";
		String	taxonomy = trans + "_tax.txt";
		int	K = 2000;
		int maxTrans = Integer.MAX_VALUE;
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
			MLTKHAlgo algo = new MLTKHAlgo();
			algo.runAlgorithm(K,dataset, taxonomy, output,true,maxTrans,true);
			algo.printStats();
		}
	}	
}
