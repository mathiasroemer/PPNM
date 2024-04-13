public class genlist<T>{
        public T[] data;
        public int size => data.Length;
        public T this[int i] => data[i];
        public genlist(){ data = new T[0]; }
        public void add(T item){
                T[] newdata = new T[size+1];
		System.Array.Copy(data,newdata,size);
                newdata[size]=item;
                data=newdata;
        }
	public static implicit operator T[](genlist<T> list){
		return list.data;
	}
	public static implicit operator genlist<T>(T[] array){
		genlist<T> list=new genlist<T>();
		list.data=array;
		return list;
	}
}
