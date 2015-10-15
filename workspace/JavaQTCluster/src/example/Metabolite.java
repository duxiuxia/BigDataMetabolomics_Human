package example;

public class Metabolite 
{
	private String Name;
	private String ID;
	private String Formula;
	private String Mass;
	private String InChIKey;
	public Metabolite(String Name, String ID, String Formula, String Mass, String InChIKey)
	{
		this.Name = Name;
		this.ID = ID;
		this.Formula = Formula;
		this.Mass = Mass;
		this.InChIKey = InChIKey;
	}
	public String getName()
	{
		return Name;
	}
	public String getID()
	{
		return ID;
	}
	public String getFormula()
	{
		return Formula;
	}
	public String getMass()
	{
		return Mass;
	}
	public String getInChIKey()
	{
		return InChIKey;
	}
}
