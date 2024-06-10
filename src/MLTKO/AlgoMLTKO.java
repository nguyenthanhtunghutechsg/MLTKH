package MLTKO;


import MLTKO_Improve.MLTKO.AlgoMLTKO_improve;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryPoolMXBean;
import java.lang.management.MemoryType;
import java.util.*;

/**
               _     _____ _____ _____ 
         _____| |___|_   _|  |  |     |
        |     | |___| | | |    -|  |  |
        |_|_|_|_|     |_| |__|__|_____|
                      top-K MLHUI Miner 

	Purpose	: This algorithm extends the TKO-Basic algorithm to mine MLHUIs
	Author	: Trinh D.D. Nguyen (dzutrinh at gmail dot com)
	Version : 1.4
	Update	: Apr, 2023
	Changes	: 
		- Apr, 04-2023: initial development
		- Apr, 11-2023: EUCP added; NUZ-item check added; code clean-up
		- Apr, 14-2023: Double utility values supported
		- Apr, 25-2023: Initial minimum utility is now properly picked (Prof. Bay suggestion).
		- Apr, 29-2023: Peak heap statistic added
						Improved initial min utility selection
						MemoryLogger removed
						Code clean-up
						Faster taxonomy access
		- Aug, 24-2023:	More comments added
**/

public class AlgoMLTKO {
	
	short version = 0x0104;					/** algorithm version number **/
	boolean useEUCP = false;				/** EUCP flag **/
	long candidateCount = 0;				/** candidate high-utility itemsets counter **/	
	double totalTime = 0; 					/** the time the algorithm terminated **/
	int huiCount = 0; 						/** the number of HUI generated  **/
	int K = 0;								/** the k parameter **/
	double minutility = 0; 					/** the internal min utility variable **/
	
	public boolean	debugging = false;		// for debugging purpose
	public boolean	topdown = true;			// top-down mining or not  
	public double 	algoRuntime = 0;		// recorded running time
	public double	algoMemUsage = 0;		// recorded memory usage 
	public boolean	benchmark = false;		// benchmark flag, if true, the algorithm will only report runtime and memory
	
	PriorityQueue<ItemsetTKO> kItemsets;	/** the top k rules found until now **/

	long transCount = 0;					/** WIP: for scalability evaluations **/
	Taxonomy taxonomy = null;				/** for describing the taxonomy of a dataset **/

	/** We create a map to store the TWU of each item */
	Map<Integer, Integer> mapItemToLevel;			// Item -> level hashmap
	Map<Integer, Double> mapItemToUtility;			// Item -> Utility
	Map<Integer, Double> mapItemToGWU;				// Item -> TWU/GWU	
	Map<Integer, List<Integer>> mapItemToAncestor;	// Real taxonomy hashmap
	Map<Integer, Map<Integer, Double>> mapFMAP;// EUCS:  key:item   key:another item   value:twu


	/** this class represent an item and its utility in a transaction */
	class Pair {		
		int item = 0;		/** an item */
		double utility = 0;	/** the utility of the item */
		public Pair(){

		}
		public Pair(int item, double utility){
			this.item = item;
			this.utility = utility;
		}/** the utility of the item */
	}

	/** 
	 * Constructor
	 */
	public AlgoMLTKO(boolean eucp) {
		this.useEUCP = eucp;
	}

	/**
	 * Run the algorithm
	 * @param input the input file path
	 * @param output the output file path
	 * @param k the parameter k
	 * @throws IOException if an error occur for reading/writing to file.
	 */
	public void runAlgorithm(String input, String tax, String output, int K)
			throws IOException {
		
		// allocate memory for the EUCP structure if requested in constructor
		if (this.useEUCP)
			mapFMAP			= new HashMap<Integer,Map<Integer, Double>>();
		
		mapItemToGWU		= new HashMap<Integer, Double>();		
		mapItemToUtility	= new HashMap<Integer, Double>();
		mapItemToLevel		= new HashMap<Integer, Integer>();
		mapItemToAncestor	= new HashMap<Integer, List<Integer>>();
		taxonomy			= new Taxonomy(tax);

		for (Map.Entry<Integer, Integer> entry : taxonomy.mapItemToParent.entrySet()) {
			Integer childItem = entry.getKey();
			Integer parentItem = entry.getValue();
			if (mapItemToAncestor.get(childItem) == null) {
				mapItemToAncestor.put(childItem, new ArrayList<>());
			}
			if (mapItemToAncestor.get(parentItem) == null) {
				mapItemToAncestor.put(parentItem, new ArrayList<>());
			}
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
		int maxLevel = 0;
		for (Map.Entry<Integer, List<Integer>> entry : mapItemToAncestor.entrySet()) {
			mapItemToLevel.put(entry.getKey(), entry.getValue().size() + 1);
			if (maxLevel < entry.getValue().size() + 1) {
				maxLevel = entry.getValue().size() + 1;
			}
		}


		long startTimestamp = System.currentTimeMillis();

		this.candidateCount = 0;
		this.minutility = 1;
		this.K = K;
		this.kItemsets = new PriorityQueue<MLTKO.ItemsetTKO>();

		// first dataset scan to calculate the TWU of each item.
		if (!benchmark) System.out.println("- First dataset scan...");

		Dataset dataset = new Dataset(input, Integer.MAX_VALUE);    // should perform similar transaction merging here, too
		transCount = dataset.getTransactions().size();

		for (int tid = 0; tid < transCount; tid++) {
			Transaction transaction = dataset.getTransactions().get(tid);
			double tu = transaction.getUtility();
			int[] items = transaction.getItems();
			double[] utils = transaction.getUtilities();
			HashMap<Integer, Double> ParentInTransaction = new HashMap<>();
			for (int i = 0; i < items.length; i++) {        // for each item, add the transaction utility to its TWU
				Integer item = items[i];
				Double util = utils[i];
				Double twu = mapItemToGWU.get(item);                        // get the current TWU of that item

				// add the utility of the item in the current transaction to its twu
				twu = (twu == null) ? tu : twu + tu;
				// calculate utility for each item on the dataset
				Double utility = mapItemToUtility.get(item);
				utility = (utility == null) ? util : utility + util;
				mapItemToUtility.put(item, utility);
				mapItemToGWU.put(item, twu);

				List<Integer> ancestors = mapItemToAncestor.get(item);
				for (Integer ancestorItem : ancestors) {
					Double utilityParent = ParentInTransaction.get(ancestorItem);
					if (utilityParent == null) {
						ParentInTransaction.put(ancestorItem, util);
					} else {
						ParentInTransaction.put(ancestorItem, utilityParent + util);
					}
				}
			}
			for (Map.Entry<Integer, Double> entry: ParentInTransaction.entrySet()) {
				Double twu = mapItemToGWU.get(entry.getKey());
				Double utility = mapItemToUtility.get(entry.getKey());
				if(twu == null) {
					mapItemToGWU.put(entry.getKey(),tu);
					mapItemToUtility.put(entry.getKey(),entry.getValue());
				}else{
					mapItemToGWU.put(entry.getKey(),twu+tu);
					mapItemToUtility.put(entry.getKey(),utility+entry.getValue());
				}
			}
		} // for tid

		// ==========================================================
		// ==== This section is based on Prof. Bay Vo suggestion ==== 
		List<Double> topUtils = new ArrayList<Double>();
		for (Integer i : mapItemToUtility.keySet()) {
			topUtils.add(mapItemToUtility.get(i));
		}
		Collections.sort(topUtils, Collections.reverseOrder());
		//System.out.println(topUtils);
		if (K >= topUtils.size())
			minutility = 0;
		else
			minutility = topUtils.get(K-1);
		if (minutility == 0) minutility = 1;
		if (!benchmark) System.out.println("  --> initial minUtil = " + minutility);
		if (!benchmark && debugging) {
			System.out.print("- Promising: ");
			for (Integer i : mapItemToUtility.keySet()) {
				System.out.print(i + " ");
			}
			System.out.println();
		}
		// ==========================================================
		
		List<List<UtilityList>> ulLists = new ArrayList<>();			// for storing ULs of items having TWU >= minutil.
		// for faster accessing the utility lists, they are are stored using map as pair: 
		// <KEY: item, VALUE: utility list associated to that item>
		Map<Integer, UtilityList> mapItemToUtilityList = new HashMap<Integer, UtilityList>();

		for (Integer item: mapItemToGWU.keySet()) {						// for each item
			if (mapItemToGWU.get(item) >= this.minutility) {			// if the item is promising (TWU >= minutil)
				UtilityList uList = new UtilityList(item);				// create an empty utility list that will be filled later.
				mapItemToUtilityList.put(item, uList);					// add the item to the list of high TWU items
			} // if
			else {
				List<Integer> listAncestorOfItem = mapItemToAncestor.get(item);
				for (int kk=0; kk < listAncestorOfItem.size(); kk++) {
					if (mapItemToGWU.get(listAncestorOfItem.get(kk)) >= this.minutility) {
						List<Integer> itemList = new ArrayList<Integer>();
						itemList.add(item);
						UtilityList tuList = new UtilityList(item);
						mapItemToUtilityList.put(item, tuList);
						break;
					} // if
				} // for k
			} // else
		} // for item

		List<List<List<Pair>>> revisedTransactions = new ArrayList<>();
		for (int i = 0; i < maxLevel; i++) {
			
			List<List<Pair>> revisedTransactionTemp = new ArrayList<>();
			List<List<Integer>> checkItemExistTemp = new ArrayList<>();
		
			for (int j = 0; j < transCount; j++) {
				List<Pair> rrTemp = new ArrayList<Pair>();
				List<Integer> ctTemp = new ArrayList<Integer>();
				revisedTransactionTemp.add(rrTemp);
				checkItemExistTemp.add(ctTemp);
			} // for j
				
			revisedTransactions.add(revisedTransactionTemp);
		} // for i
		
		if (!this.benchmark) {
			System.out.println("==== DATASET CHARACTERISTICS ====");		
			System.out.println(" Dataset: <" + input + " / " + tax + ">");			
			System.out.println(" |D|    : " + transCount);			
			System.out.println(" |GI|   : " + taxonomy.parentCount());
			System.out.println(" Depth  : " + maxLevel);			
			System.out.println(" T_max  : " + dataset.getMaxTransLength());
			System.out.println(" T_avg  : " + dataset.getAvgTransLength());
			System.out.println("=================================");
			
			System.out.println("- Second dataset scan...");
		}

		for (int tid = 0; tid < transCount; tid++) {
			Transaction transaction = dataset.getTransactions().get(tid);
    		int[] items = transaction.getItems();
    		double[] utils = transaction.getUtilities();
    		double[] remainingUtility = new double[maxLevel];
			double[] newTWU = new double[maxLevel];

			Map<Integer,Double> parentInTransaction = new HashMap<>();
			for (int i = 0; i < items.length; i++) {
				Integer item = items[i];
				if(mapItemToGWU.get(item)>=minutility){
					Pair pair = new Pair();
					pair.item = item;
					pair.utility = utils[i];
					int level = mapItemToLevel.get(item);
					revisedTransactions.get(level - 1).get(tid).add(pair);
				}
				List<Integer> listAllAncestors = mapItemToAncestor.get(item);
				for (Integer itemParent : listAllAncestors) {
					if(mapItemToGWU.get(itemParent)>=minutility){
						Double utilityParent = parentInTransaction.get(itemParent);
						if(utilityParent==null){
							parentInTransaction.put(itemParent,utils[i]);
						}else{
							parentInTransaction.put(itemParent,utilityParent+utils[i]);
						}
					}
				}
			}
			for(Map.Entry<Integer,Double> entry : parentInTransaction.entrySet()){
				int level = mapItemToLevel.get(entry.getKey())-1;
				remainingUtility[level]+=entry.getValue();
				newTWU[level]+=entry.getValue();
				revisedTransactions.get(level).get(tid).add(new Pair(entry.getKey(),entry.getValue()));
			}

			for (int i = 0; i < maxLevel; i++) {				// sort the transactions
				Collections.sort(revisedTransactions.get(i).get(tid), (o1, o2) -> compareItems(o1.item, o2.item));
			} // for i

			for(int levels = maxLevel-1; levels >= 0; levels--) {
				for(int i = 0; i < revisedTransactions.get(levels).get(tid).size(); i++) {
					Pair pair = revisedTransactions.get(levels).get(tid).get(i);

					remainingUtility[levels] = remainingUtility[levels] - pair.utility;		// subtract the utility of this item from the remaining utility
					UtilityList utilityListOfItem = mapItemToUtilityList.get(pair.item);	// get the utility list of this item

					// add new element to the utility list of this item corresponding to this transaction
					Element element = new Element(tid, pair.utility, remainingUtility[levels]);
					if (utilityListOfItem != null) utilityListOfItem.addElement(element);

					// ===============
					// EUCP
					if (this.useEUCP) {
						Map<Integer, Double> mapFMAPItem = mapFMAP.get(pair.item);
						if (mapFMAPItem == null) {
							mapFMAPItem = new HashMap<Integer, Double>();
							mapFMAP.put(pair.item, mapFMAPItem);
						} // if

						for(int j = i+1; j < revisedTransactions.get(levels).get(tid).size(); j++){
							Pair pairAfter = revisedTransactions.get(levels).get(tid).get(j);
							Double twuSum = mapFMAPItem.get(pairAfter.item);
							if(twuSum == null)
								mapFMAPItem.put(pairAfter.item, newTWU[levels]);
							else
								mapFMAPItem.put(pairAfter.item, twuSum + newTWU[levels]);
						} // for j
					}
					// ===============

				} // for i
			} // for level
		} // for tid
		
		if (!benchmark && useEUCP)
			if (debugging) {
				for (Integer i : mapItemToUtility.keySet()) {
					HashMap<Integer, Double> mapFMAPItem = (HashMap<Integer, Double>) mapFMAP.get(i);
					for (Integer j : mapItemToUtility.keySet()) {
						Double v = mapFMAPItem.get(j);
						if (v != null)
							System.out.print("["+i+","+j+"]=" + v + " ");
					}
					System.out.println();
				}
			}

		// dataset and taxonomy is now no longer needed, discard them 
		dataset = null;
		taxonomy = null;

		if (!benchmark) System.out.println("- Constructing utility lists for " + maxLevel + " level(s)...");

		for(int i = 0; i < maxLevel; i++) {
			List<UtilityList> UtilityListOfILevel = new ArrayList<>();
			for (Integer item: mapItemToGWU.keySet()){				
				if (mapItemToGWU.get(item) >= this.minutility){	// if the item is promising  (TWU >= minUtil)					
					if (mapItemToLevel.get(item) == i+1) {						
						UtilityList uList = mapItemToUtilityList.get(item);	// create an empty Utility List that we will fill later.
						UtilityListOfILevel.add(uList);	// add the item to the list of high TWU items
					} // if
				} // if
			} // for item
			
			ulLists.add(UtilityListOfILevel);
			
			// sort the list based on item's TWU in ascending order
			Collections.sort(ulLists.get(i), new Comparator<UtilityList>(){
				public int compare(UtilityList o1, UtilityList o2) {
					return compareItems(o1.item, o2.item);		// compare the TWU of the items
				}});
		} // for i

		// Mine the database recursively
		if (!benchmark) System.out.println("- MLHUI mining...");

		//
		//	The order is reversed compared to the published paper
		//          [X]		level 0
		//         /   \
		//       [a]   [b]	level 1
		//
				
		if (topdown) {
			for(int level = 0; level < maxLevel; level++) {			// top-down
								
//				if (debugging) {
//					System.out.println("-------------------");
//					System.out.println("Level = " + level + " | minUtil = " + minutility);
//					for (UtilityList ul : ulLists.get(level))
//						System.out.println(ul.toString());
//				}
	
				// mine the database recursively, levelwise
				search(new int[0], null, ulLists.get(level));
			} // for level
		}
		else {
			for(int level = maxLevel-1; level >= 0; level--) {	// bottom-up
				
//				if (debugging) {
//					System.out.println("-------------------");
//					System.out.println("Level = " + level + " | minUtil = " + minutility);
//					for (UtilityList ul : ulLists.get(level))
//						ul.print();
//				}
	
				// mine the database recursively, levelwise
				search(new int[0], null, ulLists.get(level));
			} // for level
		}
			 
		if (output != null) {
			System.out.println("- Saving results...");
			writeResultTofile(output);
		}

		totalTime = (System.currentTimeMillis() - startTimestamp) / 1000.0;

		algoRuntime = totalTime;
		algoMemUsage = peakHeapUsage();
		
		if (!benchmark) System.out.println("- Done.");
	}
	
	/**
	 * This is the recursive method to find all high utility itemsets. It writes the itemsets to the output file.
	 * 
	 * @param prefix	 This is the current prefix. Initially, it is empty.
	 * @param pUL		 This is the Utility List of the prefix. Initially, it is empty.
	 * @param ULs		 The utility lists corresponding to each extension of the prefix.
	 * @param minUtility The minUtility threshold.
	 * @throws IOException
	 */
	private void search(int[] prefix, UtilityList pUL, List<UtilityList> ULs) throws IOException {

		// For each extension X of prefix P
		for (int i = 0; i < ULs.size(); i++) {
			UtilityList X = ULs.get(i);
			
			// If pX is a high utility itemset, we save the itemset: pX
			if (X.sumIutils >= minutility)
				writeOut(prefix, X.item, X.sumIutils);

			// Z-item pruning
			if (X.sumRutils == 0) continue;
			
			// If sum remaining utilities of pX >= minUtility, we explore extensions of pX.
			// Fig 5, Line 06 from the TKO paper (this is the pruning condition)
			if (X.sumRutils + X.sumIutils >= minutility) {
				
				// This list will contain the utility lists of pX extensions.
				List<UtilityList> exULs = new ArrayList<UtilityList>();
				
				// For each extension of p appearing after X according to the ascending order
				for (int j = i + 1; j < ULs.size(); j++) {
					UtilityList Y = ULs.get(j);
					
					// ==== EUCP ====
					if (this.useEUCP) {
						Map<Integer, Double> mapTWUF = mapFMAP.get(X.item);
						if(mapTWUF != null) {
							Double twuF = mapTWUF.get(Y.item);
							if(twuF == null || twuF < minutility) continue;
						}
					}
					// ==============

					candidateCount++;
					
					// we construct the extension pXY and add it to the list of extensions of pX
					exULs.add(construct(pUL, X, Y));
				}
				
				// We create new prefix pX
				int[] newPrefix = new int[prefix.length + 1];
				System.arraycopy(prefix, 0, newPrefix, 0, prefix.length);
				newPrefix[prefix.length] = X.item;

				// We make a recursive call to discover all itemsets with the prefix pX
				search(newPrefix, X, exULs);
			}

		}
	}

	/**
	 * Method to dynamically increase minutil and maintain the list of top-K MLHUIs
	 * @param a prefix itemset
	 * @param an item to be appended to the prefix
	 * @param utility the utility of the prefix concatenated with the item
	 */
	private void writeOut(int[] prefix, int item, double utility) {
		ItemsetTKO itemset = new ItemsetTKO(prefix, item, utility);
		kItemsets.add(itemset);
		if (kItemsets.size() > this.K) {
			if (utility > this.minutility) {
				ItemsetTKO lower;
				do {
					lower = kItemsets.peek();
					if (lower == null) break; // / IMPORTANT
					kItemsets.remove(lower);
				} while (kItemsets.size() > this.K);
				this.minutility = kItemsets.peek().utility;
				//System.out.println(String.format("%.5f", this.minutility));
			}
		}
	}

	/**
	 * This method constructs the utility list of pXY
	 * @param P :  the utility list of prefix P.
	 * @param px : the utility list of pX
	 * @param py : the utility list of pY
	 * @return the utility list of pXY
	 */
	private UtilityList construct(UtilityList P, UtilityList px, UtilityList py) {
		UtilityList pxyUL = new UtilityList(py.item);	// create an empy utility list for pXY
		int xSize = px.elements.size();
		
		for (int m = 0; m < xSize; m++) {	// for each element in the utility list of pX
			Element ex = px.elements.get(m);

			// find element ey in py with tid = ex.tid
			Element ey = findElementWithTID(py, ex.tid);
			if(ey == null) continue;
						
			if(P == null){	// if the prefix p is null				
				Element eXY = new Element(	ex.tid, 	// Create the new element
											ex.iutils + ey.iutils, 
											ey.rutils);				
				pxyUL.addElement(eXY);					// add the new element to the utility list of pXY
				
			} else {
				// find the element in the utility list of p wih the same tid
				Element e = findElementWithTID(P, ex.tid);
				if(e != null){					
					Element eXY = new Element(	ex.tid,	// Create new element 
												ex.iutils + ey.iutils - e.iutils, 
												ey.rutils);					
					pxyUL.addElement(eXY);				// add the new element to the utility list of pXY
				}
			}	
		}		
		return pxyUL;									// return the utility list of pXY.
	}
	
	/**
	 * Do a binary search to find the element with a given tid in a utility list
	 * @param ulist the utility list
	 * @param tid  the tid
	 * @return  the element or null if none has the tid.
	 */
	private Element findElementWithTID(UtilityList ulist, int tid){
		List<Element> list = ulist.elements;
		
        int first = 0, last = list.size() - 1;       
        while( first <= last ) {
        	int m = ( first + last ) >>> 1;
            if(list.get(m).tid < tid) first = m + 1;
            else 
            	if(list.get(m).tid > tid) last = m - 1; 
            	else
            		return list.get(m);
        }
		return null;
	}

	/**
	 * Write the result to a file
	 * @param path the output file path
	 * @throws IOException if an exception for reading/writing to file
	 */
	public void writeResultTofile(String path) throws IOException {
		BufferedWriter writer = new BufferedWriter(new FileWriter(path));
		Iterator<ItemsetTKO> iter = kItemsets.iterator();
		while (iter.hasNext()) {
			StringBuffer buffer = new StringBuffer();
			ItemsetTKO itemset = (ItemsetTKO) iter.next();
			
			// append the prefix
			for (int i = 0; i < itemset.getItemset().length; i++) {
				buffer.append(itemset.getItemset()[i]);
				buffer.append(' ');
			}
			buffer.append(itemset.item);
			
			// append the utility value
			buffer.append(" #UTIL: ");
			buffer.append(String.format("%.5f", itemset.utility));
			
			// write to file
			writer.write(buffer.toString());
			if(iter.hasNext()){
				writer.newLine();
			}
		}
		writer.close();
	}

	private int compareItems(int item1, int item2) {
		int compare = (int) (mapItemToGWU.get(item1) - mapItemToGWU.get(item2));
		// if the same, use the lexical order otherwise use the TWU
		return (compare == 0) ? item1 - item2 : compare;
	}

	public double peakHeapUsage()	// return the peak heap usage. 
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
		} catch (Throwable t) { System.err.println("Exception in agent: " + t); }
		return retVal;
	}		
	
	/**
	 * Print statistics about the latest execution to System.out.
	 */
	public void printStats() {
		System.out.println("=============  MLTKO - " + getVersion() + " =============");
		System.out.println(" K           : " + this.K);
		System.out.println(" EUCP enabled: " + (this.useEUCP ? "yes" : "no"));
		System.out.println(" Top-down    : " + (this.topdown ? "yes" : "no"));
		System.out.println(" MLHUIs      : " + kItemsets.size());
		System.out.println(" Candidates  : " + candidateCount);
		System.out.println(" Runtime     : ~" + String.format("%.3f", totalTime) + " s");
		System.out.println(" Peak memory : ~" + String.format("%.3f", peakHeapUsage()) + " MB");
		System.out.println("==========================================");
	}
	
	/** most useless stuffs evah **/
	public String getVersion() {
		return "v." + ((version & 0xFF00) >> 8) + "." + (version & 0xFF);
	}
	
	public void printLogo() {
		System.out.println("       _     _____ _____ _____"); 
		System.out.println(" _____| |___|_   _|  |  |     |");
		System.out.println("|     | |___| | | |    -|  |  |");
		System.out.println("|_|_|_|_|     |_| |__|__|_____|");
		System.out.println("                          " + getVersion());
	}	
}
