using System;
using static System.Console;

public class genlist<T>{
	public T[] data;
	public int size = 0,capacity = 8;
	public T this[int i] => data[i];
	public genlist(){data = new T[capacity];}
	public void add(T item){
		if(size==capacity){
			T[] newdata = new T[capacity*=2];
			Array.Copy(data,newdata,size);
			data = newdata;
		}
		data[size] = item;
		size++;
	}
	public void remove(int i){
		T[] newdata = new T[size-1];
		var idx = 0;
		for(int j = 0; j < size; j++) {
			if(i==j)continue;

			newdata[idx] = data[j];
			idx++;
		}
		data = newdata;
		size--;
	}
}

public class node<T>{
	public T item;
	public node<T> next;
	public node(T item){this.item = item;}
}

public class list<T>{
        public node<T> first=null,current=null;
        public void add(T item){
                if(first==null){
                        first=new node<T>(item);
                        current=first;
                }
                else{
                        current.next = new node<T>(item);
                        current=current.next;
                }
        }
        public void start(){ current=first; }
        public void next(){ current=current.next; }
}
