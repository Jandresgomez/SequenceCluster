package main;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import ngsep.sequences.DNAShortKmer;

public class KMersTable {
	private HashMap<DNAShortKmer, Integer> table;
	private int count;
	
	public KMersTable() {
		table = new HashMap<>();
		count = 0;
	}
	
	public void add(DNAShortKmer seq) {
		Integer k = table.get(seq);
		if(k == null) {
			k = 0;
		}
		
		k++; count++;
		table.put(seq, k);
	}
	
	public void status() {
		System.out.println("Llaves: " + table.keySet().size());
		System.out.println("Count: " + count);
	}
	
	public String statusText() {
		String k = "Llaves: " + table.keySet().size() + "\n";
		k+= "Count: " + count;
		return k;
	}
	
	public DNAShortKmer getRandom() {
		Set<DNAShortKmer> s = table.keySet();
		DNAShortKmer k = (DNAShortKmer) s.toArray()[(int) (s.size()*Math.random())];
		return k;
	}
	
	public Iterator<DNAShortKmer> getKeyIterator() {
		return table.keySet().iterator();
	}
}
