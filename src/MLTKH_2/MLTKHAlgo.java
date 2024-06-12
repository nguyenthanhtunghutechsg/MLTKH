package MLTKH_2;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryType;
import java.util.*;




/* This file is copyright (c) 2012-2015 Souleymane Zida & Philippe Fournier-Viger
* 
* This file is part of the SPMF DATA MINING SOFTWARE
* (http://www.philippe-fournier-viger.com/spmf).
* 
* SPMF is free software: you can redistribute it and/or modify it under the
* terms of the GNU General Public License as published by the Free Software
* Foundation, either version 3 of the License, or (at your option) any later
* version.
* SPMF is distributed in the hope that it will be useful, but WITHOUT ANY
* WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
* A PARTICULAR PURPOSE. See the GNU General Public License for more details.
* You should have received a copy of the GNU General Public License along with
* SPMF. If not, see <http://www.gnu.org/licenses/>.
*/

/**
 * This is an implementation of the EFIM algorithm for
 * mining high-utility itemsets from a transaction database.
 * More information on the EFIM algorithm can be found in that paper: <br\>
 *
 * @author Souleymane Zida, Philippe Fournier-Viger using some code by Alan Souza
 */
public class MLTKHAlgo {

	/** the set of high-utility itemsets */
    //private Itemsets highUtilityItemsets;
    
	/** object to write the output file */
	BufferedWriter writer = null;
	
	/** the number of high-utility itemsets found (for statistics) */
	private int patternCount; 

	/** the start time and end time of the last algorithm execution */
	long startTimestamp;
	long endTimestamp;
	
	/** the minutil threshold */
	double minUtil;
	int K;
	int transCount;
	int maxLevel;
	/** if this variable is set to true, some debugging information will be shown */
    final boolean  DEBUG = false;
	
    /** The following variables are the utility-bins array 
	// Recall that each bucket correspond to an item */
    /** utility bin array for sub-tree utility */
	private double[] utilityBinArraySU;
	/** utility bin array for local utility */
	private double[] utilityBinArrayLU;
	private double[] utilityBinArrayUtility;

	/** a temporary buffer */
	private int [] temp= new int [500];
	
	/** The total time spent for performing intersections */
	long timeIntersections;
	/** The total time spent for performing database reduction */
	long timeDatabaseReduction;
	/** The total time spent for identifying promising items */
	long timeIdentifyPromisingItems;
	/** The total time spent for sorting */
	long timeSort;
	/** The total time spent for binary search */
	long timeBinarySearch;
	

	/** an array that map an old item name to its new name */
    int[] oldNameToNewNames;
    /** an array that map a new item name to its old name */
    int[] newNamesToOldNames;
    /** the number of new items */
    int newItemCount;

    /** if true, transaction merging will be performed by the algorithm */
    boolean activateTransactionMerging;

    /** A parameter for transaction merging*/
    final int MAXIMUM_SIZE_MERGING = 1000;
    
    /** number of times a transaction was read */
    long transactionReadingCount;
    /** number of merges */
    long mergeCount;

    /** number of itemsets from the search tree that were considered */
	private long candidateCount;

	/** If true, sub-tree utility pruning will be performed */
	private boolean activateSubtreeUtilityPruning;

	Map<Integer, Integer> mapItemToLevel;             // Item -> level hashmap
	double TUs[];          // Item -> Utility 	// Item -> TWU/GWU
	Map<Integer, List<Integer>> mapItemToAncestor;    // Real taxonomy hashmap
	PriorityQueue<ItemsetTKO> kItemsets; 			  //the top k rules found until now
	Taxonomy taxonomy = null;						/** for describing the taxonomy of a dataset **/
	Map<Integer, List<Integer>> mapItemToChild;
	List<Set<Integer>> listItemHUI;
	List<List<Integer>> listItemEachLevel;
	/** 
	 * Constructor
	 */
    public MLTKHAlgo() {
         
    }

    /**
     * Run the algorithm
     * @param K  the minimum utility threshold (a positive integer)
     * @param inputPath  the input file path
     * @param outputPath  the output file path to save the result or null if to be kept in memory
     * @param activateTransactionMerging 
     * @param activateSubtreeUtilityPruning 
     * @param maximumTransactionCount
       * @return the itemsets or null if the user choose to save to file
     * @throws IOException if exception while reading/writing to file
     */
    public ItemsetTKO runAlgorithm(int K, String inputPath, String taxonomyPath, String outputPath, boolean activateTransactionMerging, int maximumTransactionCount, boolean activateSubtreeUtilityPruning) throws IOException {
    	
    	// reset variables for statistics
    	mergeCount=0;
    	transactionReadingCount=0;
		timeIntersections = 0;
		timeDatabaseReduction = 0;
    	
    	// save parameters about activating or not the optimizations
    	this.activateTransactionMerging = activateTransactionMerging;
    	this.activateSubtreeUtilityPruning = activateSubtreeUtilityPruning;

		listItemHUI = new ArrayList<>();
		mapItemToLevel = new HashMap<Integer, Integer>();
		mapItemToAncestor = new HashMap<Integer, List<Integer>>();
		taxonomy = new Taxonomy(taxonomyPath);
		mapItemToChild = new HashMap<>();
		listItemEachLevel = new ArrayList<>();


		for (Map.Entry<Integer, Integer> entry : taxonomy.mapItemToParent.entrySet()) {
			Integer childItem = entry.getKey();
			Integer parentItem = entry.getValue();
			if (mapItemToAncestor.get(childItem) == null) {
				mapItemToAncestor.put(childItem, new ArrayList<>());
			}
			if (mapItemToAncestor.get(parentItem) == null) {
				mapItemToAncestor.put(parentItem, new ArrayList<>());
			}
			if (mapItemToChild.get(parentItem) == null) {
				mapItemToChild.put(parentItem, new ArrayList<>());

			}
			mapItemToChild.get(parentItem).add(childItem);
			List<Integer> ancestors = new ArrayList<>();
			Integer currentItem = childItem;
			while (currentItem != null) {
				Integer ancestor = taxonomy.mapItemToParent.get(currentItem);
				if (ancestor != null) {
					ancestors.add(ancestor);
				}
				currentItem = ancestor;
			}
			mapItemToAncestor.put(childItem, ancestors);
		}
		maxLevel = 0;
		for (Map.Entry<Integer, List<Integer>> entry : mapItemToAncestor.entrySet()) {
			mapItemToLevel.put(entry.getKey(), entry.getValue().size() + 1);
			if (maxLevel < entry.getValue().size() + 1) {
				maxLevel = entry.getValue().size() + 1;
			}
		}
		for (int i = 0; i < maxLevel; i++) {
			Set<Integer> list = new HashSet<>();
			listItemHUI.add(list);
			listItemEachLevel.add(new ArrayList<>());
			for(Map.Entry<Integer,Integer> entry  : mapItemToLevel.entrySet()){
				if(entry.getValue() - 1 == i){
					listItemEachLevel.get(i).add(entry.getKey());
				}
			}
		}
		this.candidateCount = 0;
		this.minUtil = 1;
		this.K = K;
		this.kItemsets = new PriorityQueue<>();

		// read the input file
		Dataset dataset = new Dataset(inputPath, maximumTransactionCount, mapItemToAncestor,mapItemToLevel,maxLevel );
		transCount = dataset.getTransactionInAllLevel().get(0).size();
		TUs = new double[transCount];

		// record the start time
		startTimestamp = System.currentTimeMillis();
		// reset the number of itemset found
		patternCount = 0;

		
		// if in debug mode, show the initial database in the console
		if(DEBUG)
		{
			System.out.println("===== Initial database === ");
			System.out.println(dataset.toString());
		}

        useUtilityBinArrayToCalculateLocalUtilityFirstTime(dataset);

		List<Double> listUtility = new ArrayList<>();
		for (int i = 0; i < utilityBinArrayUtility.length; i++) {
			listUtility.add(utilityBinArrayUtility[i]);
		}
		listUtility.sort(Comparator.reverseOrder());
		if(listUtility.size()>=K){
			minUtil = listUtility.get(K-1);
		}else{
			minUtil = 1;
		}

		for (int level = 0; level < maxLevel; level++) {
			List<Integer> itemsToKeep = new ArrayList<Integer>();
			for(Integer item: listItemEachLevel.get(level)){
				if(utilityBinArrayLU[item] >= minUtil) {
					itemsToKeep.add(item);
				}
			}
			insertionSort(itemsToKeep, utilityBinArrayLU);
			oldNameToNewNames = new int[dataset.getMaxItem() + 1];
			newNamesToOldNames = new int[dataset.getMaxItem() + 1];
			int currentName = 1;
			for (int j=0; j< itemsToKeep.size(); j++)
			{
				int item = itemsToKeep.get(j);
				oldNameToNewNames[item] = currentName;
				newNamesToOldNames[currentName] = item;
				itemsToKeep.set(j, currentName);
				currentName++;
			}
			newItemCount = itemsToKeep.size();
			utilityBinArraySU = new double[newItemCount + 1];
			List<Transaction> datasetInLevel = dataset.getTransactionInAllLevel().get(level);
			for(int tid=0; tid< datasetInLevel.size();tid++)
			{
				Transaction transaction  = datasetInLevel.get(tid);
				transaction.removeUnpromisingItems(oldNameToNewNames);
			}
			Collections.sort(datasetInLevel, (t1, t2) -> {
                int pos1 = t1.items.length - 1;
                int pos2 = t2.items.length - 1;
                if(t1.items.length < t2.items.length){
                    while(pos1 >=0){
                        int subtraction = t2.items[pos2]  - t1.items[pos1];
                        if(subtraction !=0){
                            return subtraction;
                        }
                        pos1--;
                        pos2--;
                    }
                    return -1;
                }else if (t1.items.length > t2.items.length){
                    while(pos2 >=0){
                        int subtraction = t2.items[pos2]  - t1.items[pos1];
                        if(subtraction !=0){
                            return subtraction;
                        }
                        pos1--;
                        pos2--;
                    }
                    return 1;
                }else{
                    while(pos2 >=0){
                        int subtraction = t2.items[pos2]  - t1.items[pos1];
                        if(subtraction !=0){
                            return subtraction;
                        }
                        pos1--;
                        pos2--;
                    }
                    return 0;
                }
            });
			int emptyTransactionCount = 0;
			for(int i=0; i< datasetInLevel.size();i++)
			{
				Transaction transaction  = datasetInLevel.get(i);
				if(transaction.items.length == 0){
					emptyTransactionCount++;
				}else{
					break;
				}
			}
			dataset.transactionInAllLevel.set(level,datasetInLevel
					.subList(emptyTransactionCount, datasetInLevel.size()));
			useUtilityBinArrayToCalculateSubtreeUtilityFirstTime(datasetInLevel);
			List<Integer> itemsToExplore = new ArrayList<Integer>();
			if(activateSubtreeUtilityPruning){
				for(Integer item : itemsToKeep){
					if (utilityBinArraySU[item] >= minUtil) {
						itemsToExplore.add(item);
					}
				}
			}
			backtrackingEFIM(datasetInLevel, itemsToKeep
					, itemsToExplore, 0);
		}

		endTimestamp = System.currentTimeMillis();

		if(outputPath != null) {
			System.out.println("- Saving results...");
			writeResultTofile(outputPath);
		}
		printPeakHeapUsage();
		return null;
    }

	/**
	 * Implementation of Insertion sort for sorting a list of items by increasing order of TWU.
	 * This has an average performance of O(n log n)
	 * @param items list of integers to be sorted
	 * @param items list the utility-bin array indicating the TWU of each item.
	 */
	public void insertionSort(List<Integer> items, double [] utilityBinArrayTWU){
		// the following lines are simply a modified an insertion sort
		for(int j=1; j< items.size(); j++){
			Integer itemJ = items.get(j);
			int i = j - 1;
			Integer itemI = items.get(i);
			double comparison;
			int levelOfItemI = mapItemToLevel.get(itemI);
			int levelOfItemJ = mapItemToLevel.get(itemJ);
			if(levelOfItemI != levelOfItemJ){
				comparison = levelOfItemI-levelOfItemJ;
			}else{
				comparison = utilityBinArrayTWU[itemI] - utilityBinArrayTWU[itemJ];
			}
			if(comparison == 0){
				comparison = itemI - itemJ;
			}
			
			while(comparison > 0){
				items.set(i+1, itemI);

				i--;
				if(i<0){
					break;
				}
				
				itemI = items.get(i);

				if(levelOfItemI != levelOfItemJ){
					comparison = levelOfItemI-levelOfItemJ;
				}else{
					comparison = utilityBinArrayTWU[itemI] - utilityBinArrayTWU[itemJ];
				}
				if(comparison == 0){
					comparison = itemI - itemJ;
				}
			}
			items.set(i+1,itemJ);
		}
	}
    
    /**
     * Recursive method to find all high-utility itemsets
     * @param transactionsOfP the list of transactions containing the current prefix P
	 * @param itemsToKeep the list of secondary items in the p-projected database
	 * @param itemsToExplore the list of primary items in the p-projected database
	 * @param prefixLength the current prefixLength
     * @throws IOException if error writing to output file
     */
    private void backtrackingEFIM( List<Transaction> transactionsOfP,
    		List<Integer> itemsToKeep, List<Integer> itemsToExplore, int prefixLength) throws IOException {
		candidateCount += itemsToExplore.size();
		for (int j = 0; j < itemsToExplore.size(); j++) {
			Integer e = itemsToExplore.get(j);
	        List<Transaction> transactionsPe = new ArrayList<>();
			int utilityPe = 0;
			Transaction previousTransaction = null;
			int consecutiveMergeCount = 0;
			long timeFirstIntersection = System.currentTimeMillis();
	        for(Transaction transaction : transactionsOfP) {
	        	long timeBinaryLocal = System.currentTimeMillis();
	        	int positionE = -1;
	        	// Variables low and high for binary search
	    		int low = transaction.offset;
	    		int high = transaction.items.length - 1;

	    		// perform binary search to find e in the transaction
	    		while (high >= low ) {
	    			int middle = (low + high) >>> 1; // divide by 2
	    			if (transaction.items[middle] < e) {
	    				low = middle + 1;
	    			}else if (transaction.items[middle] == e) {
	    				positionE =  middle;
	    				break;
	    			}  else{
	    				high = middle - 1;
	    			}
	    		}
	    		// record the time spent for performing the binary search
	        	timeBinarySearch +=  System.currentTimeMillis() - timeBinaryLocal;
	            if (positionE > -1  ) {

	            	// optimization: if the 'e' is the last one in this transaction,
	            	// we don't keep the transaction
					if(transaction.getLastPosition() == positionE){
						// but we still update the sum of the utility of P U {e}
						utilityPe  += transaction.utilities[positionE] + transaction.prefixUtility;
					}else{
						// otherwise
		            	if(activateTransactionMerging && MAXIMUM_SIZE_MERGING >= (transaction.items.length - positionE)){
			            	// we cut the transaction starting from position 'e'
							Transaction projectedTransaction = new Transaction(transaction, positionE);
							utilityPe  += projectedTransaction.prefixUtility;

							// if it is the first transaction that we read
							if(previousTransaction == null){
								// we keep the transaction in memory
								previousTransaction = projectedTransaction;
							}else if (isEqualTo(projectedTransaction, previousTransaction)){
								// If it is not the first transaction of the database and
								// if the transaction is equal to the previously read transaction,
								// we will merge the transaction with the previous one

								// increase the number of consecutive transactions merged
								mergeCount++;

								// if the first consecutive merge
								if(consecutiveMergeCount == 0){
									// copy items and their profit from the previous transaction
									int itemsCount = previousTransaction.items.length - previousTransaction.offset;
									int[] items = new int[itemsCount];
									System.arraycopy(previousTransaction.items, previousTransaction.offset, items, 0, itemsCount);
									double[] utilities = new double[itemsCount];
									System.arraycopy(previousTransaction.utilities, previousTransaction.offset, utilities, 0, itemsCount);

									// make the sum of utilities from the previous transaction
							    	int positionPrevious = 0;
									int positionProjection = projectedTransaction.offset;
									while(positionPrevious < itemsCount){
										utilities[positionPrevious] += projectedTransaction.utilities[positionProjection];
										positionPrevious++;
										positionProjection++;
									}

									// make the sum of prefix utilities
									double sumUtilities = previousTransaction.prefixUtility += projectedTransaction.prefixUtility;

									// create the new transaction replacing the two merged transactions
									previousTransaction = new Transaction(items, utilities, previousTransaction.transactionUtility + projectedTransaction.transactionUtility);
									previousTransaction.prefixUtility = sumUtilities;

								}else{
									// if not the first consecutive merge

									// add the utilities in the projected transaction to the previously
									// merged transaction
							    	int positionPrevious = 0;
									int positionProjected = projectedTransaction.offset;
									int itemsCount = previousTransaction.items.length;
									while(positionPrevious < itemsCount){
										previousTransaction.utilities[positionPrevious] += projectedTransaction.utilities[positionProjected];
										positionPrevious++;
										positionProjected++;
									}

									// make also the sum of transaction utility and prefix utility
									previousTransaction.transactionUtility += projectedTransaction.transactionUtility;
									previousTransaction.prefixUtility += projectedTransaction.prefixUtility;
								}
								// increment the number of consecutive transaction merged
								consecutiveMergeCount++;
							}else{
								// if the transaction is not equal to the preceding transaction
								// we cannot merge it so we just add it to the database
								transactionsPe.add(previousTransaction);
								// the transaction becomes the previous transaction
								previousTransaction = projectedTransaction;
								// and we reset the number of consecutive transactions merged
								consecutiveMergeCount = 0;
							}
						}else{
			            	// Otherwise, if merging has been deactivated
							// then we just create the projected transaction
							Transaction projectedTransaction = new Transaction(transaction, positionE);
							// we add the utility of Pe in that transaction to the total utility of Pe
							utilityPe  += projectedTransaction.prefixUtility;
							// we put the projected transaction in the projected database of Pe
							transactionsPe.add(projectedTransaction);
						}
					}
					// This is an optimization for binary search:
					// we remember the position of E so that for the next item, we will not search
					// before "e" in the transaction since items are visited in lexicographical order
		            transaction.offset = positionE;
	            }else{
					// This is an optimization for binary search:
					// we remember the position of E so that for the next item, we will not search
					// before "e" in the transaction since items are visited in lexicographical order
	            	transaction.offset = low;
	            }
	        }
	        // remember the total time for peforming the database projection
	        timeIntersections += (System.currentTimeMillis() - timeFirstIntersection);

	        // Add the last read transaction to the database if there is one
	        if(previousTransaction != null){
	        	transactionsPe.add(previousTransaction);
	        }

	        // Append item "e" to P to obtain P U {e}
	        // but at the same time translate from new name of "e"  to its old name
	        temp[prefixLength] = newNamesToOldNames[e];

	        // if the utility of PU{e} is enough to be a high utility itemset
	        if(utilityPe  >= minUtil)
	        {
	        	// output PU{e}
	        	output(prefixLength, utilityPe );
	        }

			//==== Next, we will calculate the Local Utility and Sub-tree utility of
	        // all items that could be appended to PU{e} ====
	        useUtilityBinArraysToCalculateUpperBounds(transactionsPe, j, itemsToKeep);

	        // we now record time for identifying promising items
			long initialTime = System.currentTimeMillis();

			// We will create the new list of secondary items
			List<Integer> newItemsToKeep = new ArrayList<Integer>();
			// We will create the new list of primary items
			List<Integer> newItemsToExplore = new ArrayList<Integer>();

			// for each item
	    	for (int k = j+1; k < itemsToKeep.size(); k++) {
	        	Integer itemk =  itemsToKeep.get(k);

	        	// if the sub-tree utility is no less than min util
	            if(utilityBinArraySU[itemk] >= minUtil) {
	            	// and if sub-tree utility pruning is activated
	            	if(activateSubtreeUtilityPruning){
	            		// consider that item as a primary item
	            		newItemsToExplore.add(itemk);
	            	}
	            	// consider that item as a secondary item
	            	newItemsToKeep.add(itemk);
	            }else if(utilityBinArrayLU[itemk] >= minUtil)
	            {
	            	// otherwise, if local utility is no less than minutil,
	            	// consider this itemt to be a secondary item
	            	newItemsToKeep.add(itemk);
	            }
	        }
	    	// update the total time  for identifying promising items
	    	timeIdentifyPromisingItems +=  (System.currentTimeMillis() -  initialTime);

			// === recursive call to explore larger itemsets
	    	if(activateSubtreeUtilityPruning){
	    		// if sub-tree utility pruning is activated, we consider primary and secondary items
	    		backtrackingEFIM(transactionsPe, newItemsToKeep, newItemsToExplore,prefixLength+1);
	    	}else{
	    		// if sub-tree utility pruning is deactivated, we consider secondary items also
	    		// as primary items
	    		backtrackingEFIM(transactionsPe, newItemsToKeep, newItemsToKeep,prefixLength+1);
	    	}
		}

    }


    /**
     * Check if two transaction are identical
     * @param t1  the first transaction
     * @param t2  the second transaction
     * @return true if they are equal
     */
    private boolean isEqualTo(Transaction t1, Transaction t2) {
    	// we first compare the transaction lenghts
		int length1 = t1.items.length - t1.offset;
		int length2 = t2.items.length - t2.offset;
		// if not same length, then transactions are not identical
    	if(length1 != length2){
    		return false;
    	}
    	// if same length, we need to compare each element position by position,
    	// to see if they are the same
    	int position1 = t1.offset;
		int position2 = t2.offset;
		
		// for each position in the first transaction
		while(position1 < t1.items.length){
			// if different from corresponding position in transaction 2
			// return false because they are not identical
			if(t1.items[position1]  != t2.items[position2]){
				return false;
			}
			// if the same, then move to next position
			position1++;
			position2++;
		}
		// if all items are identical, then return to true
		return true;
	}

	/**
	 * Scan the initial database to calculate the local utility of each item
	 * using a utility-bin array
	 * @param dataset the transaction database
	 */
	public void useUtilityBinArrayToCalculateLocalUtilityFirstTime(Dataset dataset) {

		// Initialize utility bins for all items
		utilityBinArrayLU = new double[dataset.getMaxItem() + 1];
		utilityBinArrayUtility = new double[dataset.getMaxItem() + 1];

		for (int level = 0; level < maxLevel; level++) {
			List<Transaction> tranInLevel = dataset.getTransactionInAllLevel().get(level);
			for (int tid = 0;tid<transCount;tid++) {
				// for each item
				Transaction transaction = tranInLevel.get(tid);
				TUs[tid] = transaction.transactionUtility;
				for (int i = 0; i < transaction.getItems().length; i++) {
					int item = transaction.getItems()[i];
					utilityBinArrayLU[item] += transaction.transactionUtility;
					utilityBinArrayUtility[item]+=transaction.getUtilities()[i];

				}
			}
		}
	}
	
	
	/**
	 * Scan the initial database to calculate the sub-tree utility of each item
	 * using a utility-bin array
	 * @param dataset the transaction database
	 */
	public void useUtilityBinArrayToCalculateSubtreeUtilityFirstTime(List<Transaction> dataset) {

		int sumSU;
		for (Transaction transaction : dataset) {
			// We will scan the transaction backward. Thus,
			// the current sub-tree utility in that transaction is zero
			// for the last item of the transaction.
			sumSU = 0;

			// For each item when reading the transaction backward
			for(int i = transaction.getItems().length-1; i >=0; i--) {
				// get the item
				Integer item = transaction.getItems()[i];

				// we add the utility of the current item to its sub-tree utility
				sumSU += transaction.getUtilities()[i];
				// we add the current sub-tree utility to the utility-bin of the item
				utilityBinArraySU[item] += sumSU;
			}
		}

	}

    /**
     * Utilize the utility-bin arrays to calculate the sub-tree utility and local utility of all
     * items that can extend itemset P U {e}
     * @param transactionsPe transactions the projected database for P U {e}
     * @param j the position of j in the list of promising items
     * @param itemsToKeep the list of promising items
     */
    private void useUtilityBinArraysToCalculateUpperBounds(List<Transaction> transactionsPe,
    		int j, List<Integer> itemsToKeep) {

    	// we will record the time used by this method for statistics purpose
		long initialTime = System.currentTimeMillis();

		// For each promising item > e according to the total order
		for (int i = j + 1; i < itemsToKeep.size(); i++) {
			Integer item = itemsToKeep.get(i);
			// We reset the utility bins of that item for computing the sub-tree utility and
			// local utility
			utilityBinArraySU[item] = 0;
			utilityBinArrayLU[item] = 0;
		}

		int sumRemainingUtility;
		// for each transaction
		for (Transaction transaction : transactionsPe) {
			// count the number of transactions read
			transactionReadingCount++;

			// We reset the sum of reamining utility to 0;
			sumRemainingUtility = 0;
			// we set high to the last promising item for doing the binary search
			int high = itemsToKeep.size() - 1;

			// for each item in the transaction that is greater than i when reading the transaction backward
			// Note: >= is correct here. It should not be >.
			for (int i = transaction.getItems().length - 1; i >= transaction.offset; i--) {
				// get the item
				int item = transaction.getItems()[i];

				// We will check if this item is promising using a binary search over promising items.

				// This variable will be used as a flag to indicate that we found the item or not using the binary search
				boolean contains = false;
				// we set "low" for the binary search to the first promising item position
				int low = 0;

				// do the binary search
				while (high >= low) {
					int middle = (low + high) >>> 1; // divide by 2
					int itemMiddle = itemsToKeep.get(middle);
					if (itemMiddle == item) {
						// if we found the item, then we stop
						contains = true;
						break;
					} else if (itemMiddle < item) {
						low = middle + 1;
					} else {
						high = middle - 1;
					}
				}
				// if the item is promising
				if (contains) {
					// We add the utility of this item to the sum of remaining utility
					sumRemainingUtility += transaction.getUtilities()[i];
					// We update the sub-tree utility of that item in its utility-bin
					utilityBinArraySU[item] += sumRemainingUtility + transaction.prefixUtility;
					// We update the local utility of that item in its utility-bin
					utilityBinArrayLU[item] += transaction.transactionUtility + transaction.prefixUtility;
				}
			}
		}
		// we update the time for database reduction for statistics purpose
		timeDatabaseReduction += (System.currentTimeMillis() - initialTime);
    }



    /**
     * Save a high-utility itemset to file or memory depending on what the user chose.
     * @param tempPosition of the temp
	 * @param utility of the itemset
     * @throws IOException if error while writting to output file
     */
    private void output(int tempPosition, int utility) throws IOException {
		int [] prefix = new int[tempPosition];
		int item = temp[tempPosition];
		for (int i = 0; i < tempPosition; i++) {
			prefix[i]=temp[i];
		}
		ItemsetTKO itemset = new ItemsetTKO(prefix, item, utility);
		kItemsets.add(itemset);
		if (utility > this.minUtil) {
			int level = mapItemToLevel.get(item) - 1;
			for (int itemInItemset : prefix) {
				listItemHUI.get(level).add(itemInItemset);
			}
			listItemHUI.get(level).add(item);
		}
		if (kItemsets.size() > this.K) {
			if (utility > this.minUtil) {
				ItemsetTKO lower;
				do {
					lower = kItemsets.peek();
					if (lower == null) break; // / IMPORTANT
					kItemsets.remove(lower);
				} while (kItemsets.size() > this.K);
				this.minUtil = kItemsets.peek().utility;
				//System.out.println(String.format("%.5f", this.minutility));
			}
		}
    }


    public static void printPeakHeapUsage()
    {
    	try {
            List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
            // we print the result in the console
			double total = 0;
			for (MemoryPoolMXBean memoryPoolMXBean : pools) {
				if (memoryPoolMXBean.getType() == MemoryType.HEAP) {
					long peakUsed = memoryPoolMXBean.getPeakUsage().getUsed();
					//System.out.println(String.format("Peak used for: %s is %.2f", memoryPoolMXBean.getName(), (double)peakUsed/1024/1024));
					total = total + peakUsed;
				}
			}
			System.out.println(String.format("Total heap peak used: %f MB", total/1024/1024));

       } catch (Throwable t) {
            System.err.println("Exception in agent: " + t);
       }
    }


    /**
     * Print statistics about the latest execution of the EFIM algorithm.
     */
	public void printStats() {

		System.out.println("========== EFIM v97 - STATS ============");
		System.out.println(" minUtil = " + minUtil);
		System.out.println(" High utility itemsets count: " + patternCount);
		System.out.println(" Total time ~: " + (endTimestamp - startTimestamp)
				+ " ms");
		System.out.println(" Transaction merge count ~: " + mergeCount);
		System.out.println(" Transaction read count ~: " + transactionReadingCount);

		// if in debug mode, we show more information
		if(DEBUG) {

			System.out.println(" Time intersections ~: " + timeIntersections
					+ " ms");
			System.out.println(" Time database reduction ~: " + timeDatabaseReduction
					+ " ms");
			System.out.println(" Time promising items ~: " + timeIdentifyPromisingItems
					+ " ms");
			System.out.println(" Time binary search ~: " + timeBinarySearch
					+ " ms");
			System.out.println(" Time sort ~: " + timeSort	+ " ms");
		}
		System.out.println(" Max memory:" +  peakHeapUsage());
		System.out.println(" Candidate count : "             + candidateCount);
		System.out.println("=====================================");
	}
	public double peakHeapUsage()    // return the peak heap usage.
	{
		double retVal = 0;
		try {
			List<MemoryPoolMXBean> pools = ManagementFactory.getMemoryPoolMXBeans();
			double total = 0;
			for (MemoryPoolMXBean memoryPoolMXBean : pools) {
				if (memoryPoolMXBean.getType() == MemoryType.HEAP) {
					long peakUsed = memoryPoolMXBean.getPeakUsage().getUsed();
					total = total + peakUsed;
				}
			}
			retVal = total / (1024L * 1024L);
		} catch (Throwable t) {
			System.err.println("Exception in agent: " + t);
		}
		return retVal;
	}

	public void writeResultTofile(String path) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(path));
		Iterator<ItemsetTKO> iter = kItemsets.iterator();
		while (iter.hasNext()) {
			StringBuffer buffer = new StringBuffer();
			ItemsetTKO itemset = (ItemsetTKO) iter.next();
			for (int i = 0; i < itemset.getItemset().length; i++) {
				buffer.append(itemset.getItemset()[i]);
				buffer.append(' ');
			}
			buffer.append(itemset.item);
			buffer.append(" #UTIL: ");
			buffer.append(String.format("%.5f", itemset.utility));

			// write to file
			writer.write(buffer.toString());
			if (iter.hasNext()) {
				writer.newLine();
			}
		}
		writer.close();
	}
}
