package MLTKH;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
//import java.util.Set;

// mlTKOs
// Class represents a Transaction Dataset 
//
// Coded by Trinh D.D. Nguyen
// Version 1.0 - Mar, 2023
//
//This class represents a quantitative database
public class Dataset {
	
	List<List<Transaction>> transactionInAllLevel;						// the list of transactions in this dataset
	private int maxItem = 0;							// the largest item name
	private long maxTransLength = 0;					// longest transaction length
	private long sumTransLength = 0;					// total transactions length

	// main constructor
    public Dataset(String datasetPath, int maximumTransactionCount,
                   Map<Integer,List<Integer>> mapItemToAncestor,
                   Map<Integer,Integer> mapItemToLevel, int maxLevel) throws IOException {
        transactionInAllLevel = new ArrayList<>();
        for (int i = 0; i < maxLevel; i++) {
            transactionInAllLevel.add(new ArrayList<Transaction>());
        }
        BufferedReader br = new BufferedReader(new FileReader(datasetPath));
        String line;
        int i = 0;

        while((line = br.readLine()) != null) { 
			if (line.isEmpty() == true || 				// bypass comments and empty lines
				line.charAt(0) == '#' ||  
				line.charAt(0) == '@') continue;
			i++;
            List<Transaction> transactionInEachLevel = createTransaction(line,mapItemToAncestor,mapItemToLevel,maxLevel);
            for (int j = 0; j < maxLevel; j++) {
                transactionInAllLevel.get(j).add(transactionInEachLevel.get(j));
            }
	    	if(i == maximumTransactionCount) break;		// prevent exceeding the number of transactions
        }
        br.close();
    }

    // create a transaction object from a string read from the input file
    private List<Transaction> createTransaction(String line,  Map<Integer,List<Integer>> mapItemToAncestor,
                                                Map<Integer,Integer> mapItemToLevel, int maxLevel) {

        String[] split = line.split(":");								// split the line into tokens using ":"
    	double transactionUtility = Double.parseDouble(split[1]);			// Get the transaction utility
        String[] itemsString = split[0].split(" ");						// Get the list of items
        String[] itemsUtilitiesString = split[2].split(" ");			// Get the list of item utilities
        int[] items = new  int[itemsString.length];						// store the items and their utilities
        double[] utilities = new double[itemsString.length];

        Map<Integer,Double> mapItemInTransaction = new HashMap<>();
        for (int i = 0; i < items.length; i++) {						// for each item        	
        	items[i] = Integer.parseInt(itemsString[i]);				// store that item        	
        	utilities[i] = Double.parseDouble(itemsUtilitiesString[i]);	// and its utility in that transaction
            List<Integer> listOfAncestors = mapItemToAncestor.get(items[i]);
            mapItemInTransaction.put(items[i],utilities[i]);
            for (Integer parent: listOfAncestors) {
                Double utility = mapItemInTransaction.get(parent);
                if( utility == null){
                    mapItemInTransaction.put(parent,utilities[i]);
                }else{
                    mapItemInTransaction.put(parent,utilities[i]+utility);
                }
            }
        }
        List<List<Integer>> itemsInEachLevel = new ArrayList<>();
        List<List<Double>> utilitiesInEachLevel = new ArrayList<>();
        for (int i = 0; i < maxLevel; i++) {
            itemsInEachLevel.add(new ArrayList<>());
            utilitiesInEachLevel.add(new ArrayList<>());
        }
        for(Map.Entry<Integer,Double> entry: mapItemInTransaction.entrySet()){
            int level = mapItemToLevel.get(entry.getKey())-1;
            itemsInEachLevel.get(level).add(entry.getKey());
            utilitiesInEachLevel.get(level).add(entry.getValue());
            if(entry.getKey() > maxItem) maxItem = entry.getKey();
        }
        List<Transaction> result = new ArrayList<>();
        for (int i = 0; i < maxLevel; i++) {
            int[] itemArray = itemsInEachLevel.get(i).stream().mapToInt(Integer::intValue).toArray();
            double[] doubleArray = utilitiesInEachLevel.get(i).stream().mapToDouble(Double::doubleValue).toArray();
            result.add(new Transaction(itemArray
                    ,doubleArray
                    , transactionUtility));
        }
        if (maxTransLength < items.length)	maxTransLength = items.length;
        sumTransLength += items.length;
		
        return result;	// create the transaction
    }

    // returns the list of all transactions
    public List<List<Transaction>> getTransactions() {
        return transactionInAllLevel;
    }

    // returns the largest item within the database
    public int getMaxItem() {
        return maxItem;
    }

    // returns the maximum transaction length
    public long getMaxTransLength() {
    	return maxTransLength;
    }

    // returns the average transaction length
    public double getAvgTransLength() {
    	return (double) sumTransLength / transactionInAllLevel.get(0).size();
    }
}
