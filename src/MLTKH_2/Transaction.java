package MLTKH_2;
 
// ML-HUI-MINER ALGORITHM MULTI-CORE
// ---------------------------------
// This algorithm extends the FHM algorithm to mine multi-level high-utility itemsets
//
// Coded by Trinh D.D. Nguyen
// Version 1.1 - Nov, 2020
//

// This class represents a transaction from the quantitative database
public class Transaction {

    public static int[] tempItems = new int[2000];/** a buffer to store items of an itemset*/
    public static double[] tempUtilities = new double[2000];/** a buffer to store utilities of an itemset */


    int[] items;				// list of items contained within the transaction
    double[] utilities;			// list of utilities associated to items within the transaction
    double transactionUtility; 	// the transaction utility value (TU)
    double prefixUtility;/** a buffer to store utilities of an itemset */
    /** an offset pointer, used by projected transactions*/
    int offset;

    // main constructor
    public Transaction(int[] items, double[] utilities, double transactionUtility) {
    	this.items = items;
    	this.utilities = utilities;
    	this.transactionUtility = transactionUtility;
        this.offset = 0;
        this.prefixUtility = 0;
    }
    public int getLastPosition(){
        return items.length -1;
    }
    /**
     * Constructor for a projected transaction
     * @param transaction the transaction that will be projected (it may be an original transaction
     * or a previously projected transaction
     * @param offsetE an offset over the original transaction for projecting the transaction
     */
    public Transaction(Transaction transaction, int offsetE) {
        // copy items and utilities from the original transaction
        this.items = transaction.getItems();
        this.utilities = transaction.getUtilities();

        // copy the utility of element e
        double utilityE = this.utilities[offsetE];

        // add the  utility of item e to the utility of the whole prefix used to project the transaction
        this.prefixUtility = transaction.prefixUtility + utilityE;

        // we will now calculate the remaining utility.
        // It is the transaction utility minus the profit of the element that was removed
        this.transactionUtility = transaction.transactionUtility - utilityE;
        // and we also need to subtract the utility of all items before e
        // but after the previous offset
        for(int i = transaction.offset; i < offsetE; i++){
            this.transactionUtility -= transaction.utilities[i];
        }
        // remember the offset for this projected transaction
        this.offset = offsetE+1;
    }

    /**
     * This method removes unpromising items from the transaction and at the same time rename
     * items from old names to new names
     * @param oldNamesToNewNames An array indicating for each old name, the corresponding new name.
     */
    public void removeUnpromisingItems(int[] oldNamesToNewNames) {
        // In this method, we used buffers for temporary storing items and their utilities
        // (tempItems and tempUtilities)
        // This is for memory optimization.

        // for each item
        int i = 0;
        for(int j=0; j< items.length;j++) {
            int item = items[j];

            // Convert from old name to new name
            int newName = oldNamesToNewNames[item];

            // if the item is promising (it has a new name)
            if(newName != 0) {
                // copy the item and its utility
                tempItems[i] = newName;
                tempUtilities[i] = utilities[j];
                i++;
            }else{
                // else subtract the utility of the item
                transactionUtility -= utilities[j];
            }
        }
        // copy the buffer of items back into the original array
        this.items = new int[i];
        System.arraycopy(tempItems, 0, this.items, 0, i);

        // copy the buffer of utilities back into the original array
        this.utilities = new double[i];
        System.arraycopy(tempUtilities, 0, this.utilities, 0, i);

        // Sort by increasing TWU values
        insertionSort(this.items, this.utilities);
    }

    /**
     * Implementation of Insertion sort for integers.
     * This has an average performance of O(n log n)
     * @param items array of integers
     */
    public static void insertionSort(int [] items,  double[] utitilies){
        for(int j=1; j< items.length; j++){
            int itemJ = items[j];
            double utilityJ = utitilies[j];
            int i = j - 1;
            for(; i>=0 && (items[i]  > itemJ); i--){
                items[i+1] = items[i];
                utitilies[i+1] = utitilies[i];
            }
            items[i+1] = itemJ;
            utitilies[i+1] = utilityJ;
        }
    }

    
    // returns the list of items in the transaction
    public int[] getItems() {
        return items;
    }
    
    // returns the list of utilities of items in the transaction
    public double[] getUtilities() {
        return utilities;
    }

    // returns the transaction utility value
	public double getUtility()
	{
		return transactionUtility;
	}

	// return the length of the transaction
	public long length()
	{
		return items.length;
	}
}
