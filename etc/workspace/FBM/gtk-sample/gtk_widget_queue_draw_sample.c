#include <gtk/gtk.h>
#include <stdlib.h>

char *imagefile = "image.jpg";

struct rectangle {
  double x0;
  double y0;
  double x1;
  double y1;
} rect = {
  50.0,
  50.0,
  80.0,
  80.0
};


static void
cb_button_left (GtkWidget *button, gpointer user_data)
{
  rect.x0 -= 10;
  rect.x1 -= 10;
  gtk_widget_queue_draw(GTK_WIDGET(user_data));
}

static void
cb_button_right (GtkWidget *button, gpointer user_data)
{
  rect.x0 += 10;
  rect.x1 += 10;
  gtk_widget_queue_draw(GTK_WIDGET(user_data));
}

static void
cb_button_top (GtkWidget *button, gpointer user_data)
{
  rect.y0 -= 10;
  rect.y1 -= 10;
  gtk_widget_queue_draw(GTK_WIDGET(user_data));
  printf("(%f,%f)\n",rect.y0, rect.y1);
}

static void
cb_button_down (GtkWidget *button, gpointer user_data)
{
  rect.y0 += 10;
  rect.y1 += 10;
  gtk_widget_queue_draw(GTK_WIDGET(user_data));
  printf("(%f,%f)\n",rect.y0, rect.y1);
}


gboolean cb_expose_event(GtkWidget *widget,
			 GdkEventExpose *event,
			 gpointer user_data)
{
  GdkWindow *drawable = widget->window;
  cairo_t *cr;
  cr = gdk_cairo_create (drawable);

  gdk_cairo_set_source_pixbuf (cr, (GdkPixbuf *)user_data, 0.0, 0.0);
  cairo_paint (cr);
  

  cairo_set_source_rgb (cr, 1.0, 0.0, 0.0);
  cairo_rectangle (cr, rect.x0, rect.y0, rect.x1, rect.y1);
  cairo_stroke (cr);

  cairo_destroy (cr);

  return FALSE;
}


int main (int argc, char *argv[])
{
  GtkWidget *window;
  GtkWidget *canvas;
  GdkPixbuf *pixbuf;

  
  gtk_init (&argc, &argv);
  pixbuf = gdk_pixbuf_new_from_file (imagefile,NULL);

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title (GTK_WINDOW (window),
			"cairo_source_image Sample");
  gtk_widget_set_size_request (window,
			       gdk_pixbuf_get_width (pixbuf),
			       gdk_pixbuf_get_height (pixbuf)*1.8);

  g_signal_connect (G_OBJECT (window), "destroy",
		    G_CALLBACK (gtk_main_quit), NULL);

  canvas = gtk_drawing_area_new ();
  g_signal_connect (G_OBJECT (canvas), "expose-event",
		    G_CALLBACK (cb_expose_event), pixbuf);
  {
    GtkWidget *vbox = gtk_vbox_new(FALSE, 0);
    {
      gtk_container_add (GTK_CONTAINER (window), vbox);
      gtk_box_pack_start(GTK_BOX(vbox), canvas, TRUE, TRUE, 0);
      {
	GtkWidget *hbox = gtk_hbox_new(FALSE, 0);
	GtkWidget *button;


	gtk_container_add (GTK_CONTAINER(vbox), hbox);

	button = gtk_button_new_with_label("left");
	g_signal_connect (G_OBJECT(button), "clicked",
			  G_CALLBACK(cb_button_left), canvas);
	gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);	
	
	button = gtk_button_new_with_label("right");
	g_signal_connect (G_OBJECT(button), "clicked",
			  G_CALLBACK(cb_button_right), canvas);
	gtk_box_pack_start(GTK_BOX(hbox), button, TRUE, TRUE, 0);	
      }
    }
  }

  gtk_widget_show_all (window);
  gtk_main ();
  g_object_unref (pixbuf);

  return 0;
}
