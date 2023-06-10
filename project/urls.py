from django.contrib import admin
from django.urls import path, include
from main import views


admin.site.site_header = "Панель управления"

urlpatterns = [
    path('admin/', admin.site.urls),
    path('main/', include('main.urls')),
    path('', views.LoginView.as_view(), name='login'),
    ]

