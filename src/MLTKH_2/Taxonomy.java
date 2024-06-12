package MLTKH_2;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

// Taxonomy
// --------
// Class to store a dataset's taxonomy on memory for quick access
// 
// Coded by Trinh D.D. Nguyen, 
// Version 1.1 - Sep 2022
public class Taxonomy {

	public class Tuple {				// representing a tuple consists of parent node and its child
		int parent = 0;
		int child = 0;
		
		Tuple(int p, int c) {			// default constructor
			parent = p;
			child = c;
		}
	}	

	public Map<Integer, Integer> mapItemToParent;
	
	// default constructor
	public Taxonomy() { 				 
		mapItemToParent = new LinkedHashMap<Integer, Integer>();
	}

	// another constructor
	public Taxonomy(String filename) throws IOException { 
		mapItemToParent = new LinkedHashMap<Integer, Integer>();
		load(filename);
	}
	
	// add a tuple to the taxonomy 
	public void add(int p, int c) {
		mapItemToParent.put(c, p);	
	}
	
	// load taxonomy from text file
	public void load(String filename) throws IOException {
		BufferedReader	reader = new BufferedReader(new FileReader(filename)); 
		String			line;

		try {
			while ((line = reader.readLine()) != null) {		// scanning through the text file
			
				if (line.isEmpty() == true || line.charAt(0)=='#' || line.charAt(0)=='@') 
					continue;									// skipping comments and empty lines
											
				String	tokens[] = line.split(",");				// splitting string using ','														
				int	child = Integer.parseInt(tokens[0]);		// child comes first								
				int parent = Integer.parseInt(tokens[1]);		// then its parent							
			
				add(parent, child);								// then add this tuple into the list
			}
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		finally {
			if(reader != null) reader.close(); 
		}
	}
	
	/*
	// access index-th parent
	public int parent(int index) {
		return taxonomy.get(index).parent;
	}
	
	// access index-th child
	public int child(int index) {
		return taxonomy.get(index).child;
	}

	// access index-th tuple
	public Tuple get(int index) {
		return taxonomy.get(index);
	}

	// return the total tuples in the taxonomy 
	public int size() {
		return taxonomy.size();
	}
	*/
	// return the number of parent nodes in the taxonomy - for statistical purposes only
	public int parentCount() {
		List<Integer> parents = new ArrayList<Integer>();
		for (Integer v : mapItemToParent.values()) {
			if (!parents.contains(v))
				parents.add(v);
		}
		return parents.size();
	}
}
